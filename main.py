#!/usr/bin/python3

import fem
from load import load_grid
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib
import os
import imageio
from matplotlib import animation
import sys

np.set_printoptions(suppress=True)
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(linewidth=300)
np.set_printoptions(precision=3)


max_H = 0.1
max_B = 0.15
min_H = 0.0
min_B = 0.0
integration_points_number = 3
n_H = 15
n_B = 15
#grid = fem.Grid(max_H, max_B, n_H, n_B, integration_points_number)

# load from file
grid, options = load_grid('grid.txt')
d_tau = options["SimulationStepTime"]
tau_end = options["SimulationTime"]


n_H = grid.elements[0][3]
n_B = n_H

max_B = float('-inf')
min_B = float('inf')
max_H = float('-inf')
min_H = float('inf')
for node in grid.nodes:
    if node.x > max_B:
        max_B = node.x
    if node.x < min_B:
        min_B = node.x
    if node.y > max_H:
        max_H = node.y
    if node.y < min_H:
        min_H = node.y


d_tau = 1.0
tau = 0.0
tau_end = 20.0

tmpe = []

y, x = np.meshgrid(np.linspace(min_H, max_H, n_H), np.linspace(min_B, max_B, n_B))

for i in np.arange(0.0, tau_end, d_tau):
    t0 = np.array([node.t for node in grid.nodes])
    H = grid.H_aggregated + grid.C_aggregated/d_tau
    P = grid.P_aggregated + np.dot(grid.C_aggregated/d_tau, t0)

    t1 = np.linalg.solve(H, P)
    
    for idx in range(len(grid.nodes)):
        grid.nodes[idx].t = t1[idx]

    tmpe  = t1

a = []
for node in grid.nodes:
    a.append(node.t)


a = np.array(a).reshape(n_B, n_H)

fig, ax = plt.subplots()
c = ax.pcolormesh(x, y, a, cmap='RdBu')

fig.colorbar(c, ax=ax)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

