#!/usr/bin/python3

import fem
import numpy as np
import math

def gauss(func, points, table):
    result = 0.0
    for key, value in table[points].items():
        result = result + func(key) * value
    return result


def gauss2d(func, points, table):
    result = 0.0
    for first_key in table[points]:
        for second_key in table[points]:
            result = result + func(first_key, second_key) * table[points][first_key] * table[points][second_key]
    return result

            
def func1(x):
    return 5*x**2 + 3*x + 6


def func2(x, y):
    return 5*x**2*y**2 + 3*x*y + 6


table = [
        {
            -1/math.sqrt(3): 1,
            1/math.sqrt(3): 1
        },
        {
            -math.sqrt(3/5): 5/9,
            0: 8/9,
            math.sqrt(3/5): 5/9
        },
        {
            -0.861136: 0.347855,
            -0.339981: 0.652145,
            0.339981: 0.652145,
            0.861136: 0.347855
        }
]


class Elem4_2D:
    def __init__(self, eta = 0, ksi = 0):
        self.ksi = ksi
        self.eta = eta



elem4_2d_eta_table = []
elem4_2d_ksi_table = []

class elem4_2d:
    def __init__(self, ksi=0, eta=0):
        self.ksi = ksi
        self.eta = eta

    def __repr__(self):
        return f'[{self.ksi}, {self.eta}]'

    def __str__(self):
        return f'[{self.ksi}, {self.eta}]'


eta_shape_funcs = [
    lambda x: -1/4 * (1-x),
    lambda x: -1/4 * (1+x),
    lambda x: 1/4 * (1+x),
    lambda x: 1/4 * (1-x),
]

ksi_shape_funcs = [
    lambda x: -1/4 * (1-x),
    lambda x: 1/4 * (1-x),
    lambda x: 1/4 * (1+x),
    lambda x: -1/4 * (1+x)
]

elems = []

schema_2_pts = []
for i in range(4):
    elem4_2d_eta_table.append([None] * 4)
    elem4_2d_ksi_table.append([None] * 4)
    elems.append([Elem4_2D()] * 4)

    keys = list(table[0].keys())
    elem4_2d_eta_table[i][0] = eta_shape_funcs[i](keys[0])
    elem4_2d_eta_table[i][3] = eta_shape_funcs[i](keys[0])
    elems[i][0].eta = eta_shape_funcs[i](keys[0])
    elems[i][3].eta = eta_shape_funcs[i](keys[0])

    elem4_2d_ksi_table[i][0] = ksi_shape_funcs[i](keys[0])
    elem4_2d_ksi_table[i][3] = ksi_shape_funcs[i](keys[0])
    elems[i][0].ksi = ksi_shape_funcs[i](keys[0])
    elems[i][3].ksi = ksi_shape_funcs[i](keys[0])

    elem4_2d_eta_table[i][1] = eta_shape_funcs[i](keys[1])
    elem4_2d_eta_table[i][2] = eta_shape_funcs[i](keys[1])
    elems[i][1].eta = eta_shape_funcs[i](keys[1])
    elems[i][2].eta = eta_shape_funcs[i](keys[1])

    elem4_2d_ksi_table[i][1] = ksi_shape_funcs[i](keys[1])
    elem4_2d_ksi_table[i][2] = ksi_shape_funcs[i](keys[1])
    elems[i][1].ksi = ksi_shape_funcs[i](keys[1])
    elems[i][2].ksi = ksi_shape_funcs[i](keys[1])


grid = fem.Grid(0.2, 0.1, 5, 4)

nodes = [(0.0,0.0), (0.025, 0.0), (0.025, 0.025), (0.0, 0.025)]
jacobian = []

elem4_2d_ksi_table = np.transpose(np.array(elem4_2d_ksi_table))
elem4_2d_eta_table = np.transpose(np.array(elem4_2d_eta_table))
print(elem4_2d_ksi_table)
print(elem4_2d_eta_table)


for e in elems:
    for r in e:
        print(r.eta, end=" ")
    print()

result_x = 0
result_y = 0
for i in range(4):
    result_x += elem4_2d_ksi_table[0][i] * nodes[i][0]
    result_y += elem4_2d_eta_table[0][i] * nodes[i][1]


j = np.array([[result_x, 0], [0, result_y]])
print(j)
print("det", np.linalg.det(j))
print("1/det", 1/np.linalg.det(j))
print(np.linalg.inv(j))

print(grid.elements)

print(grid.nodes)

for i in range(grid.nE):
    j = 0
    for j in range(grid.elements[i].nodes):
        pass
        


def jacobian(i, j, jacobian, jacobian_inv, elem2d_4, grid):
    pass


























































"""
elem9_2d_eta_table = []

for i in range(4):
    elem9_2d_eta_table.append([None] * 9)
    keys = list(table[1].keys())
    elem9_2d_eta_table[i][0] = eta_shape_funcs[i](keys[0])
    elem9_2d_eta_table[i][3] = eta_shape_funcs[i](keys[0])
    elem9_2d_eta_table[i][6] = eta_shape_funcs[i](keys[0])

    elem9_2d_eta_table[i][1] = eta_shape_funcs[i](keys[1])
    elem9_2d_eta_table[i][4] = eta_shape_funcs[i](keys[1])
    elem9_2d_eta_table[i][7] = eta_shape_funcs[i](keys[1])

    elem9_2d_eta_table[i][2] = eta_shape_funcs[i](keys[2])
    elem9_2d_eta_table[i][5] = eta_shape_funcs[i](keys[2])
    elem9_2d_eta_table[i][8] = eta_shape_funcs[i](keys[2])

print()
print("eta 9 elem")
for i in range(len(elem9_2d_eta_table)):
    print(elem9_2d_eta_table[i])

print()
print("ksi 9 elem")
"""
