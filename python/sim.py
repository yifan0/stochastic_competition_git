import numpy as np
from numpy.random import binomial
import random
import matplotlib.pyplot as plt

size = 100
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
    speciation_event_count = binomial(world_size, specrate)
    speciation_events = random.sample(range(world_size), speciation_event_count)
    for index in speciation_events:
        i = int(index / size)
        j = int(index % size)
        ratio = 1 + random.random() * mutsize
        if random.random() < 0.5:
            ratio = 1.0 / ratio
        psuccess = p * ratio / (p * (ratio - 1) + 1)
        if random.random() <= psuccess:
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
                    
    invasion_updates = {}
    # Do invasion events
    for i in range(size):
        for j in range(size):
            if land_mask[i, j] > 0:
                neighbors = []
                for x in range(-1,2,1):
                    row = (i + x) % size
                    for y in range(-1,2,1):
                        col = (j + y) % size
                        neighbors.append(land_grid[row, col])
                del neighbors[4] # remove self
                next_val = random.choices(neighbors, weights=neighbors, k=1)[0]
                invasion_updates[(i,j)] = next_val

    for pair in invasion_updates:
        i,j = pair
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

plt.pcolor(land_grid, cmap=plt.cm.seismic, vmin=land_grid.min(), vmax=land_grid.max())
plt.colorbar()
plt.savefig(f"heatmap.png")
