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

shape_funcs = [
        lambda ksi, eta: 0.25*(1-ksi)*(1-eta),
        lambda ksi, eta: 0.25*(1+ksi)*(1-eta),
        lambda ksi, eta: 0.25*(1+ksi)*(1+eta),
        lambda ksi, eta: 0.25*(1-ksi)*(1+eta)
]

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

points = [
        [-1/math.sqrt(3), -1/math.sqrt(3)],
        [1/math.sqrt(3), -1/math.sqrt(3)],
        [1/math.sqrt(3), 1/math.sqrt(3)],
        [-1/math.sqrt(3), 1/math.sqrt(3)],
        ]

elems = []

schema_2_pts = []
for i in range(4):
    elem4_2d_eta_table.append([None] * 4)
    elem4_2d_ksi_table.append([None] * 4)
    elems.append([])
    for j in range(4):
        elems[i].append(Elem4_2D())

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

elems = np.transpose(np.array(elems))
for i in range(len(elems[2])):
    elems[2][i].ksi, elems[3][i].ksi = elems[3][i].ksi, elems[2][i].ksi
    elems[1][i].eta, elems[3][i].eta = elems[3][i].eta, elems[1][i].eta


grid = fem.Grid(0.2, 0.1, 30, 20, 2)

nodes = [(0.0,0.0), (0.025, 0.0), (0.025, 0.025), (0.0, 0.025)]
jacobian = []

elem4_2d_ksi_table = np.transpose(np.array(elem4_2d_ksi_table))
elem4_2d_eta_table = np.transpose(np.array(elem4_2d_eta_table))



result_x = 0
result_y = 0
for i in range(4):
    result_x += elem4_2d_ksi_table[0][i] * nodes[i][0]
    result_y += elem4_2d_eta_table[0][i] * nodes[i][1]


j = np.array([[result_x, 0], [0, result_y]])


jacobians = []

"""
for element in grid.elements:
    jacobians.append([])
    for i in range(len(element)):
        pcx = np.array([x.eta for x in elems[i]])
        pcy = np.array([x.ksi for x in elems[i]])
        x = grid.getXCoords(element)
        y = grid.getYCoords(element)
        result_x = np.sum(pcx*x)
        result_y = np.sum(pcy*y)
        result_xx = np.sum(pcy*x)
        result_yy = np.sum(pcx*y)
        jacobian = np.array([[result_xx,result_x], [result_y, result_yy]])
        jacobians[-1].append(jacobian)

jacobians = np.array(jacobians)
"""


"""
for i in range(len(jacobians)):
    for j in range(len(jacobians[i])):
        e = elems[j]
        for k in range(len(e)):
            t = np.transpose(np.matrix([e[k].ksi, e[k].eta]))
            inv = np.linalg.inv(jacobians[i][j])
            result = np.dot(inv, t)
            grid.elements[i].x_derivatives[j][k] = result[0][0]
            grid.elements[i].y_derivatives[j][k] = result[1][0]

arr = np.zeros((4,4))
for i in range(4):
    arr[i] = shape_values = np.array([
            shape_funcs[0](*points[i]),
            shape_funcs[1](*points[i]),
            shape_funcs[2](*points[i]),
            shape_funcs[3](*points[i])
            ])


for i in range(len(grid.elements)):
    C = np.zeros((4,4))
    for j in range(len(jacobians[i])):
        Nx = np.dot(grid.elements[i].x_derivatives[j,:].reshape(4,1), grid.elements[i].x_derivatives[j,:].reshape(1, -1))
        Ny = np.dot(grid.elements[i].y_derivatives[j,:].reshape(4,1), grid.elements[i].y_derivatives[j,:].reshape(1, -1))
        H = 25 * (Nx + Ny) * np.linalg.det(jacobians[i][j])
        grid.elements[i].H[j] = H
        C = np.zeros((4,4))
        element = grid.elements[i]
        element.C_sum += grid.c * grid.ro * np.dot(arr[j].reshape(4,1), arr[j].reshape(1, -1)) * np.linalg.det(jacobians[i][j])


sides_coords = np.array([
    [[-1/math.sqrt(3), -1], [1/math.sqrt(3), -1]],
    [[1, -1/math.sqrt(3)], [1, 1/math.sqrt(3)]],
    [[1/math.sqrt(3), 1], [-1/math.sqrt(3), 1]],
    [[-1, 1/math.sqrt(3)], [-1, -1/math.sqrt(3)]]
    ])

detJ = np.array([
    (grid.B/(grid.nB-1))/2,
    (grid.H/(grid.nH-1))/2
    ])

for element in grid.elements:
    H = np.zeros((4,4))
    

for element in grid.elements:
    for i in range(len(element.sides)):
        if not np.array_equal(element.sides[i], np.array([0,0])):
            values_1 = np.zeros((4))
            values_2 = np.zeros((4))
            values_1[i%4] = shape_funcs[i%4](*sides_coords[i%4][0])
            values_2[i%4] = shape_funcs[i%4](*sides_coords[i%4][1])

            values_1[(i+1)%4] = shape_funcs[(i+1)%4](*sides_coords[i%4][0])
            values_2[(i+1)%4] = shape_funcs[(i+1)%4](*sides_coords[i%4][1])

            point_1 = np.dot(values_1.reshape(4,1), values_1.reshape(1,-1))
            point_2 = np.dot(values_2.reshape(4,1), values_2.reshape(1,-1))
            H = 300 * (point_1 + point_2) * detJ[i%2]
            element.H_BC[i] = H


# calculate H aggregated

for element in grid.elements:
    element.sum_H()
    H = element.H_sum
    ids_matrix = np.zeros((4,4,2))
    for i in range(len(element.nodes)):
        for j in range(len(element.nodes)):
            ids_matrix[i][j][0] = element.nodes[i]
            ids_matrix[j][i][1] = element.nodes[i]
    for i in range(len(ids_matrix)):
        for j in range(len(ids_matrix[i])):
            H_BC_sum = np.zeros((4,4))
            for H_BC in element.H_BC:
                H_BC_sum += H_BC
            grid.H_aggregated[int(ids_matrix[i][j][0])][int(ids_matrix[i][j][1])] += element.H_sum[i][j] + H_BC_sum[i][j]
            grid.C_aggregated[int(ids_matrix[i][j][0])][int(ids_matrix[i][j][1])] += element.C_sum[i][j]



# calculate P vector v2

ambient_temperature = 1200

for element in grid.elements:
    for i in range(len(element.sides)):
        if not np.array_equal(element.sides[i], np.array([0,0])):
            values_1 = np.zeros((4))
            values_2 = np.zeros((4))
            values_1[i%4] = shape_funcs[i%4](*sides_coords[i%4][0])
            values_2[i%4] = shape_funcs[i%4](*sides_coords[i%4][1])

            values_1[(i+1)%4] = shape_funcs[(i+1)%4](*sides_coords[i%4][0])
            values_2[(i+1)%4] = shape_funcs[(i+1)%4](*sides_coords[i%4][1])

            element.P += (300 * (values_1.reshape(4,1)*ambient_temperature + (values_2.reshape(4,1)*ambient_temperature)) * detJ[i%2]).reshape(4)


for element in grid.elements:
    P_local = element.P
    for i in range(len(element.nodes)):
        grid.P_aggregated[element.nodes[i]] += element.P[i]

"""

d_tau = 5.0
tau = 0.0
tau_end = 300.0


tmpe = []

for i in np.arange(0.0, tau_end, d_tau):
    t0 = np.array([node.t for node in grid.nodes])
    H = grid.H_aggregated + grid.C_aggregated/d_tau
    P = grid.P_aggregated + np.dot(grid.C_aggregated/d_tau, t0)

    t1 = np.linalg.solve(H, P)

    
    for idx in range(len(grid.nodes)):
        grid.nodes[idx].t = t1[idx]

    tmpe = np.flip(t1.reshape((grid.nB, grid.nH))).T
    


plot = plt.matshow(tmpe, cmap='OrRd')
plt.colorbar()

plt.show()
