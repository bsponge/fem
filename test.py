import math
import numpy as np

def cartesian_product(x):
    result = set()
    for first_value in x:
        for second_value in x:
            result.add((first_value, second_value))
    return result


gauss_values  = [
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
gauss_points = [list(cartesian_product(gauss_values[i].keys()))
        for i in range(len(gauss_values))]

i = 0


def sides_sort(x):
    x.sort(key=lambda x: x[1])
    length = len(gauss_values[i])
    cnt = 0
    b = True
    while cnt < len(x):
        if b:
            slc = x[cnt:cnt+length]
            slc.sort(key=lambda x: x[0])
            x[cnt:cnt+length] = slc
        else:
            slc = x[cnt:cnt+length]
            slc.sort(key=lambda x: -x[0])
            x[cnt:cnt+length] = slc
        cnt += length
        b = not b
    for v in x:
        print(v)


l = list(gauss_points[i])
sides_sort(l)
