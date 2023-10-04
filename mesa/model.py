import mesa
from numpy.random import binomial
import random

class Plant(mesa.Agent):
    def __init__(self, unique_id, model):
        super().__init__(unique_id, model)
        self.lifehistory = 1
        self.next_val = 1
    
    def speciate(self, ratio):
        self.lifehistory = self.lifehistory * ratio

    def step(self):
        neighbors = self.model.grid.get_neighbors(self.pos, moore=True, include_center=False)
        self.next_val = random.choices(neighbors, weights=[agent.lifehistory for agent in neighbors], k=1)[0].lifehistory

    def advance(self):
        self.lifehistory = self.next_val

class PlantModel(mesa.Model):
    def __init__(self, width, height, individuals_per_patch = 10, specrate = 0.0001, invrate = 0.2, mutsize = 0.1):
        self.num_agents = width*height
        self.specrate = specrate
        self.invrate = invrate
        self.mutsize = mutsize
        self.grid = mesa.space.SingleGrid(width, height, True)
        self.schedule = mesa.time.SimultaneousActivation(self)
        for x in range(width):
            for y in range(height):
                a = Plant(y*width + x, self)
                self.schedule.add(a)
                self.grid.place_agent(a, (x, y))

        self.p = 1 / individuals_per_patch
        self.timescale = 100*(self.num_agents*1.0*self.p)

    def step(self):
        # Speciation events
        speciation_events = binomial(self.num_agents, self.specrate)
        speciation_agents = random.sample(range(self.num_agents), speciation_events)
        for agent_id in speciation_agents:
            ratio = 1 + random.random() * self.mutsize
            if random.random() < 0.5:
                ratio = 1.0 / ratio
            psuccess = self.p * ratio / (self.p*(ratio-1)+1)
            if random.random() <= psuccess:
                self.schedule._agents[agent_id].speciate(ratio)

        # Call agents to do invasion events
        self.schedule.step()
