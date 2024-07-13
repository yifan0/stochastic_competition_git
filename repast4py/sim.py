from typing import Dict, Tuple
from mpi4py import MPI
from repast4py import core, schedule, random, space, parameters
from repast4py.space import DiscretePoint as dpt
from repast4py import context as ctx
import numpy as np
import matplotlib.pyplot as plt
import math

# Usage: mpirun -n 4 python sim.py sim.yaml

size = 0
specrate = 0
invrate = 0
mutsize = 0
p = 0
output_file = "test.out"
OFFSETS = np.array([-1, 1])
rng : np.random.Generator = np.random.default_rng()

class Plant(core.Agent):
    TYPE = 0
    lifehistory = 1
    
    def __init__(self, local_id: int, rank: int, pt: dpt):
        super().__init__(id=local_id, type=Plant.TYPE, rank=rank)
        self.pt = pt
        self.lifehistory = 1
    
    def save(self) -> Tuple:
        return (self.uid, self.lifehistory)
        #return (self.uid, self.pt, self.lifehistory)
    
    def speciation(self):
        # Has been selected for speciation from the binomial distribution
        ratio = 1 + rng.random() * mutsize
        if rng.random() < 0.5:
            ratio = 1.0 / ratio
        psuccess = p * ratio / (p * (ratio - 1) + 1)
        if rng.random() <= psuccess:
            lifehistory *= ratio
    
    def invasion(self, context : ctx):
        neighbors = []
        for i in range(-1,2,1):
            x = (self.pt.x + i) % size
            for j in range(-1,2,1):
                y = (self.pt.y + j) % size
                if not (i == 0 and j == 0):
                    agent_id = (x*size+y, self.type, self.uid_rank) # NOTE: the id is Tuple(id, type, rank)
                    neighbors.append(context.agent(agent_id).lifehistory)
        neighbors_dist = np.asarray(neighbors)
        neighbors_dist = neighbors_dist / np.sum(neighbors_dist)
        next_val = rng.choice(neighbors, p=neighbors_dist)
        self.lifehistory = next_val

plant_cache = {}

def restore_plant(plant_data: Tuple):
    uid = plant_data[0]
    if uid in plant_cache: # TODO: how do we know this is up to date with lifehistory value?
        plant_cache[uid].lifehistory = plant_data[1]
        return plant_cache[uid]
    else:
        pt = space.DiscretePoint(math.floor(uid[0]/size), uid[0]%size)
        plant = Plant(uid[0], uid[2], pt)
        plant.lifehistory = plant_data[1]
        plant_cache[uid] = plant
        return plant

class Model:
    rank = -1
    def __init__(self, comm: MPI.Intracomm, params: Dict):
        global rng
        rng = np.random.default_rng(params['random.seed'])
        # create the schedule
        self.runner = schedule.init_schedule_runner(comm)
        self.runner.schedule_repeating_event(1, 1, self.speciation)
        self.runner.schedule_repeating_event(1, 1, self.invasion)
        timescale = int(100 * (size*1.0 / params['p']))
        self.runner.schedule_stop(timescale)
        self.runner.schedule_end_event(self.at_end)

        # create the context
        self.context = ctx.SharedContext(comm)

        # create a bounding box equal to the size of the entire global world grid
        box = space.BoundingBox(0, size, 0, size, 0, 0)
        # create a SharedGrid of 'box' size with sticky borders that allows multiple agents
        # in each grid location.
        self.grid = space.SharedGrid(name='grid', bounds=box, borders=space.BorderType.Periodic,
                                     occupancy=space.OccupancyType.Single, buffer_size=2, comm=comm)
        self.context.add_projection(self.grid)

        rank = comm.Get_rank()
        for i in range(size):
            for j in range(size):
                pt = space.DiscretePoint(i, j)
                plant = Plant(j*size+i, rank, pt)
                self.context.add(plant)
                self.grid.move(plant, pt)


        # TODO: add logging -- maybe for the phylogenetic tree?

    def speciation(self):
        speciation_event_count = rng.binomial(size*size, specrate)
        for i in range(speciation_event_count):
            self.context.agents()[i].speciation(self.grid) 
        #self.context.synchronize(restore_plant) # TODO: correct and reactivate
        # tick = self.runner.schedule.tick
        # self.data_set.log(tick) # Don't need to log every event that happens

    def invasion(self):
        # TODO: if we need a mask it'll need to be a property of the Plant object
        for plant in self.context.agents():
            plant.invasion(self.context)
        # tick = self.runner.schedule.tick
        # self.data_set.log(tick)

    def at_end(self):
        # print each agent's value
        # TODO: how to do this in distributed mem?
        with open(output_file + str(self.rank) + ".csv", 'w') as fout:
            for i in range(size):
                for j in range(size):
                    pt = space.DiscretePoint(i, j)
                    fout.write(str(self.grid.get_agent(pt)))
                    if j < size-1:
                        fout.write(", ")
                if i < size-1:
                    fout.write("\n")
        # TODO: if we need a mask it'll need to be a property of the Plant object
        # TODO: close anything that needs to be closed
        pass

    def start(self):
        self.runner.execute()

def run(params: Dict):
    global size, specrate, invrate, mutsize, p, output_file
    size = params['world.size']
    specrate = params['specrate']
    invrate = params['invrate']
    mutsize = params['mutsize']
    p = params['p']
    output_file = params['output_file']
    model = Model(MPI.COMM_WORLD, params)
    model.start()


if __name__ == "__main__":
    parser = parameters.create_args_parser()
    args = parser.parse_args()
    params = parameters.init_params(args.parameters_file, args.parameters)
    run(params)
