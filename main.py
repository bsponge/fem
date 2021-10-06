import fem
import numpy as np

grid = fem.Grid(0.5, 0.1, 5, 4)
print(grid.nodes)
print(grid.elements)
