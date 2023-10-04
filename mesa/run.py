#!python

import matplotlib.pyplot as plt
from model import PlantModel
import numpy as np
import seaborn as sns

model = PlantModel(100, 100)
for i in range(100000):
    model.step()

# analysis
agent_vals = np.zeros((model.grid.width, model.grid.height))
for cell, (x,y) in model.grid.coord_iter():
    agent_vals[x][y] = cell.lifehistory

g = sns.heatmap(agent_vals, cmap="rocket", annot=True, cbar=False, square=True)
plt.savefig("mesa.png")
