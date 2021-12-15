#!/usr/bin/python3

import fem
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
np.set_printoptions(linewidth=350)
np.set_printoptions(precision=3)

def myplot(x, y, s, bins=1000):
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=bins)
    heatmap = gaussian_filter(heatmap, sigma=s)

    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    return heatmap.T, extent


height = 0.1
width = 0.1
integration_points_number = 2
n_H = 31
n_B = 31
grid = fem.Grid(height, width, n_H, n_B, integration_points_number)


d_tau = 1.0
tau = 0.0
tau_end = 20.0


tmpe = []

size_x = n_B + (n_B-1)*(integration_points_number-2)
size_y = n_H + (n_H-1)*(integration_points_number-2)


y, x = np.meshgrid(np.linspace(0.0, width, size_x), np.linspace(0.0, height, size_y))

for i in np.arange(0.0, tau_end, d_tau):
    t0 = np.array([node.t for node in grid.nodes])
    H = grid.H_aggregated + grid.C_aggregated/d_tau
    P = grid.P_aggregated + np.dot(grid.C_aggregated/d_tau, t0)

    t1 = np.linalg.solve(H, P)

    
    for idx in range(len(grid.nodes)):
        grid.nodes[idx].t = t1[idx]

    tmpe  = t1

"""
print(x.shape)
a = np.zeros(x.shape)

fig, ax = plt.subplots()
c = ax.pcolormesh(x, y, a, cmap='RdBu')

fig.colorbar(c, ax=ax)
plt.show()
"""
print(t1)
print(max(t1))
