import numpy as np
from numpy import random
import matplotlib.pyplot as plt
from numba import jit
from numba import types
from numba.typed import Dict

# NOTE: does not allow wrap around, assumes boundary values are 0

# Fix to bypass np.random.choice probability option
@jit
def rand_choice_nb(arr, prob):
    """
    :param arr: A 1D numpy array of values to sample from.
    :param prob: A 1D numpy array of probabilities for the given samples.
    :return: A random sample from the given array with a given probability.
    """
    return arr[np.searchsorted(np.cumsum(prob), np.random.random(), side="right")]

@stencil
def update_mask(land):
    local_max = max(land[-1,-1], land[-1,0], land[-1,1], land[0,-1], land[0,1], land[1,-1], land[1,0], land[1,1])
    return local_max * p / (local_max * p + land[0,0] * (1-p))

@jit
def sim():
    size : int = 20
    world_size = size*size
    specrate = 0.0001
    mutsize = 0.1
    p = 0.1
    timescale = int(100 * (size * 1.0 / p))
    invrate = 0.2

    land_grid = np.ones((size, size))
    land_mask = np.zeros((size, size))

    for step in range(timescale):
        # Invasion events
        speciation_event_count = np.random.binomial(world_size, specrate)
        speciation_events = np.random.choice(world_size, size=speciation_event_count, replace=False)
        # TODO: reframe as vectorized
        for index in speciation_events:
            i = int(index / size)
            j = int(index % size)
            ratio = 1 + np.random.random() * mutsize
            if np.random.random() < 0.5:
                ratio = 1.0 / ratio
            psuccess = p * ratio / (p * (ratio - 1) + 1)
            if np.random.random() <= psuccess:
                land_grid[i,j] *= ratio
                # update the mask
                for x in range(-1,2,1):
                    row = (i + x) % size
                    for y in range(-1,2,1):
                        col = (j + y) % size
                        local_max = land_grid[row, col]
                        for xx in range(-1,2,1):
                            neighbor_row = (row + xx) % size
                            for yy in range(-1,2,1):
                                neighbor_col = (col + yy) % size
                                local_max = max(local_max, land_grid[neighbor_row, neighbor_col])
                        land_mask[row, col] = local_max * p / (local_max * p + land_grid[row,col] * (1-p))
                        
        invasion_updates = Dict.empty(key_type=types.int64, value_type=types.float64)
        #invasion_updates = {}
        # Do invasion events
        # TODO: redo as vectorized
        for i in range(size):
            for j in range(size):
                if land_mask[i, j] > 0:
                    neighbors = np.zeros((3,3))
                    for x in range(-1,2,1):
                        row = (i + x) % size
                        for y in range(-1,2,1):
                            col = (j + y) % size
                            neighbors[x+1, y+1] = land_grid[row, col]
                    neighbors[1,1] = 0 # may not select self
                    neighbor_weights = neighbors / neighbors.sum() # normalize to unit vector
                    next_val = rand_choice_nb(neighbors.flatten(), neighbor_weights.flatten())
                    invasion_updates[i*size+j] = next_val

        for pair in invasion_updates:
            i = int(pair / size)
            j = pair % size
            land_mask[i,j] = invasion_updates[pair]
            for x in range(-1,2,1):
                row = (i + x) % size
                for y in range(-1,2,1):
                    col = (j + y) % size
                    local_max = 0
                    for xx in range(-1,2,1):
                        row = (i + x + xx) % size
                        for yy in range(-1,2,1):
                            col = (j + y + yy) % size
                            local_max = max(local_max, land_grid[neighbor_row, neighbor_col])
                    land_mask[row, col] = local_max * p / (local_max * p + land_grid[row,col] * (1-p))
                    row = (i + x) % size
                    col = (j + y) % size
                    if local_max == 0:
                        land_mask[row, col] = 0
                    else:
                        land_mask[row, col] = local_max * p / (local_max * p + land_grid[row,col] * (1-p))

    #print(land_grid)
    np.savetxt("sim.csv", land_grid, delimiter=",")
    #plt.pcolor(land_grid, cmap=plt.cm.seismic, vmin=land_grid.min(), vmax=land_grid.max())
    #plt.colorbar()
    #plt.savefig(f"heatmap.png")

sim()
