#!/usr/bin/python

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
    def __init__(self, eta, ksi):
        self.ksi = ksi
        self.eta = eta


elem4_2d_table = []

shape_funcs = [
    lambda x: -1/4 * (1-x),
    lambda x: -1/4 * (1+x),
    lambda x: 1/4 * (1+x),
    lambda x: 1/4 * (1-x),
]


schema_2_pts = []
for i in range(4):
    elem4_2d_table.append([None] * 4)
    keys = list(table[0].keys())
    elem4_2d_table[i][0] = shape_funcs[i](keys[0])
    elem4_2d_table[i][3] = shape_funcs[i](keys[0])

    elem4_2d_table[i][1] = shape_funcs[i](keys[1])
    elem4_2d_table[i][2] = shape_funcs[i](keys[1])

for i in range(len(elem4_2d_table)):
    print(elem4_2d_table[i])


elem9_2d_table = []

for i in range(4):
    elem4_2d_table.append([None] * 4)
    keys = list(table[0].keys())
    elem4_2d_table[i][0] = shape_funcs[i](keys[0])
    elem4_2d_table[i][3] = shape_funcs[i](keys[0])

    elem4_2d_table[i][1] = shape_funcs[i](keys[1])
    elem4_2d_table[i][2] = shape_funcs[i](keys[1])
