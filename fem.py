from typing import List
import numpy as np
import math


class Node:
    def __init__(self, x: float = 0.0, y: float = 0.0):
        self.x = x
        self.y = y
        self.BC = False
        self.t = 100.0

    def __repr__(self):
        return f'[{self.x:.3f}, {self.y:.3f}]'

class Side:
    def __init__(self, points: np.ndarray, nodes: List, idx: int, horizontal: bool):
        self.points = points
        self.nodes = nodes
        self.horizontal = horizontal
        self.idx = idx


class Element:
    integration_points_number = 2

    def __init__(self, nodes: np.ndarray = np.zeros(4, dtype=int)):
        self.nodes = nodes
        self.x_derivatives = np.zeros((4,4))
        self.y_derivatives = np.zeros((4,4))
        self.H = np.zeros((Element.integration_points_number**2, len(self.nodes), len(self.nodes)))
        self.H_sum = np.zeros((len(self.nodes), len(self.nodes)))
        self.H_BC = np.zeros((4, len(self.nodes), len(self.nodes)))
        self.P = np.zeros((len(self.nodes)))
        self.C = np.zeros((len(self.nodes), len(self.nodes)))
        self.C_sum = np.zeros((len(self.nodes), len(self.nodes)))
        self.dN_dx = np.zeros((Element.integration_points_number**2, len(self.nodes)))
        self.dN_dy = np.zeros((Element.integration_points_number**2, len(self.nodes)))
        self.sides = []


    def __repr__(self):
        return self.nodes.__repr__()

    def __getitem__(self, index):
        return self.nodes[index]

    def __len__(self):
        return len(self.nodes)

    def sum_H(self):
        for H in self.H:
            self.H_sum += H


class Grid:

    @staticmethod
    def empty():
        grid = Grid()
        grid.integration_points_number = 3
        gauss_points = [list(cartesian_product(Grid.gauss_values[i].keys()))
                for i in range(len(Grid.gauss_values))]
        for gauss_pts in gauss_points:
            sides_sort(gauss_pts, Grid.gauss_values, grid.integration_points_number-2)

        grid.gauss_2d_points = gauss_points[grid.integration_points_number-2]
        grid.shape_functions = Grid.shape_funcs
        grid.integration_points_number = grid.integration_points_number

        grid.deta_derivatives_values = np.zeros(
                (len(grid.gauss_2d_points),
                    len(Grid.shape_funcs))
                )
        grid.dksi_derivatives_values = np.zeros(
                (len(grid.gauss_2d_points),
                    len(Grid.shape_funcs))
                )

        grid.elements = []
        grid.nodes = []
        grid.ambient_temperature = 1200
        
        # calculate derivatives
        for i in range(len(Grid.shape_funcs_derivatives[0])):
            for j in range(len(grid.gauss_2d_points)):
                grid.deta_derivatives_values[j][i] = Grid.shape_funcs_derivatives[0][i](*grid.gauss_2d_points[j])
                grid.dksi_derivatives_values[j][i] = Grid.shape_funcs_derivatives[1][i](*grid.gauss_2d_points[j])

        return grid

    def __init__(self, h: float = 0.0,
            b: float = 0.0,
            n_h: float = 0,
            n_b: float = 0,
            integration_points_number: int = 3):
        Element.integration_points_number = integration_points_number

        if h == 0 or b == 0 or n_h == 0 or n_b == 0:
            return

        gauss_points = [list(cartesian_product(Grid.gauss_values[i].keys()))
                for i in range(len(Grid.gauss_values))]
        for gauss_pts in gauss_points:
            sides_sort(gauss_pts, Grid.gauss_values, integration_points_number-2)

        self.gauss_2d_points = gauss_points[integration_points_number-2]
        self.shape_functions = Grid.shape_funcs
        self.integration_points_number = integration_points_number

        self.deta_derivatives_values = np.zeros(
                (len(self.gauss_2d_points),
                    len(Grid.shape_funcs))
                )
        self.dksi_derivatives_values = np.zeros(
                (len(self.gauss_2d_points),
                    len(Grid.shape_funcs))
                )


        # calculate derivatives
        for i in range(len(Grid.shape_funcs_derivatives[0])):
            for j in range(len(self.gauss_2d_points)):
                self.deta_derivatives_values[j][i] = Grid.shape_funcs_derivatives[0][i](*self.gauss_2d_points[j])
                self.dksi_derivatives_values[j][i] = Grid.shape_funcs_derivatives[1][i](*self.gauss_2d_points[j])


        self.c = 700
        self.ro = 7800
        self.alpha = 300.0
        self.k = 25
        self.ambient_temperature = 1200.0
        self.H = h
        self.B = b
        self.nH = n_h
        self.nB = n_b
        self.nN = self.nH * self.nB
        self.nE = (self.nH - 1) * (self.nB - 1)
        x = self.nB
        self.elements = []
        self.nodes = []

        dx = self.B / (self.nB-1)
        dy = self.H / (self.nH-1)

        # create nodes
        for x in range(self.nB):
            for y in range(self.nH):
                x_cord = x*dx
                y_cord = y*dy
                if x == self.nB-1:
                    x_cord = self.B
                if y == self.nH-1:
                    y_cord = self.H
                
                
                node = Node(x_cord, y_cord)
                if x_cord == 0.0 or x_cord == self.B:
                    node.BC = True
                if y_cord == 0.0 or y_cord == self.H:
                    node.BC = True
                self.nodes.append(node)


        self.nodes = np.array(self.nodes)
       

        # init indexes
        for i in range(0, self.nB-1):
            for j in range(0, self.nH-1):
                idx = []
                # lower side
                idx.append(i*self.nH+j)
                idx.append((i+1)*self.nH + j)
                idx.append((i+1)*self.nH + j + 1)
                idx.append(i*self.nH+j+1)
                

                element = Element(idx)
                self.elements.append(element)


        self.init_aggregated_matrices()

        self.calculate_jacobians()
        self.calculate_dNs()
        
        
        self.calculate_values()
        
        self.calculate_H_and_C()
        self.load_sides()

        self.calculate_H_BC()

        self.calculate_P()
        self.aggregate()

    def init_aggregated_matrices(self):
        self.H_aggregated = np.zeros((len(self.nodes), len(self.nodes)))
        self.P_aggregated = np.zeros((len(self.nodes)))
        self.C_aggregated = np.zeros((len(self.nodes), len(self.nodes)))



    def calculate_dNs(self):
        for i in range(len(self.jacobians)):
            for j in range(len(self.jacobians[i])):
                dksi = self.dksi_derivatives_values[j]
                deta = self.deta_derivatives_values[j]
                for k in range(len(dksi)):
                    t = np.transpose(np.matrix([
                        dksi[k],
                        deta[k]
                        ]))
                    inv = np.linalg.inv(self.jacobians[i][j])
                    result = np.dot(inv, t)
                    self.elements[i].dN_dx[j][k] = result[0][0]
                    self.elements[i].dN_dy[j][k] = result[1][0]




    def calculate_values(self):
        self.values_arr = np.zeros((len(self.gauss_2d_points), len(self.elements[0].nodes)))
        for i in range(len(self.values_arr)):
            for j in range(len(self.shape_functions)):
                self.values_arr[i][j] = self.shape_functions[j](*self.gauss_2d_points[i])



    def calculate_H_and_C(self):
        for i in range(len(self.elements)):
            for j in range(len(self.jacobians[i])):
                Nx = np.dot(self.elements[i].dN_dx[j,:].reshape(self.elements[i].dN_dx.shape[1], 1),
                        self.elements[i].dN_dx[j,:].reshape(1, -1))
                Ny = np.dot(self.elements[i].dN_dy[j,:].reshape(self.elements[i].dN_dy[j,:].shape[0], 1),
                        self.elements[i].dN_dy[j,:].reshape(1, -1))
                H = self.k * (Nx + Ny) * np.linalg.det(self.jacobians[i][j]) * Grid.gauss_values[self.integration_points_number-2][self.gauss_2d_points[j][0]] * Grid.gauss_values[self.integration_points_number-2][self.gauss_2d_points[j][1]]
                self.elements[i].H[j] = H
                element = self.elements[i]
                element.C_sum += self.c * self.ro * Grid.gauss_values[self.integration_points_number-2][self.gauss_2d_points[j][0]] * Grid.gauss_values[self.integration_points_number-2][self.gauss_2d_points[j][1]] * np.dot(self.values_arr[j].reshape(len(self.values_arr[j]), 1), self.values_arr[j].reshape(1, -1)) * np.linalg.det(self.jacobians[i][j])

            for j in range(len(self.jacobians[i])):
                self.elements[i].H_sum += self.elements[i].H[j]



    
    def calculate_H_BC(self):
        for element in self.elements:
            for side in element.sides:
                values = np.zeros((self.integration_points_number, len(element.nodes)))
                for i in range(self.integration_points_number):
                    for j in range(len(self.shape_funcs)):
                        values[i][j] = self.shape_funcs[j](*side.points[i])

                axis = 'x'
                if not side.horizontal:
                    axis = 'y'

                minimum_node = self.min_node(side.nodes, axis)
                maximum_node = self.max_node(side.nodes, axis)
                det_J = pythagoras(minimum_node, maximum_node)/2

                coords = []
                for point in side.points:
                    if point[0] != 1.0 and point[0] != -1.0:
                        coords.append(point[0])
                    if point[1] != 1.0 and point[1] != -1.0:
                        coords.append(point[1])


                for i in range(len(values)):
                    element.H_BC[side.idx] += self.alpha * (Grid.gauss_values[self.integration_points_number-2][coords[i]]*np.dot(values[i].reshape(len(element.nodes), 1), values[i].reshape(1, -1))) * det_J


            

    def calculate_P(self):
        for element in self.elements:
            for side in element.sides:
                values = np.zeros((self.integration_points_number, len(element.nodes)))
                for i in range(self.integration_points_number):
                    for j in range(len(self.shape_funcs)):
                        values[i][j] = self.shape_funcs[j](*side.points[i])


                axis = 'x'
                if not side.horizontal:
                    axis = 'y'

                minimum_node = self.min_node(side.nodes, axis)
                maximum_node = self.max_node(side.nodes, axis)
                det_J = pythagoras(minimum_node, maximum_node)/2

                coords = []
                for point in side.points:
                    if point[0] != 1.0 and point[0] != -1.0:
                        coords.append(point[0])
                    if point[1] != 1.0 and point[1] != -1.0:
                        coords.append(point[1])

                for i in range(len(values)):
                    element.P += (self.alpha * (Grid.gauss_values[self.integration_points_number-2][coords[i]]*values[i].reshape(len(element.nodes), 1)*self.ambient_temperature) * det_J).reshape(len(element.nodes))




    def aggregate(self):
        for element in self.elements:
            P_local = element.P
            for i in range(len(element.nodes)):
                self.P_aggregated[element.nodes[i]] += element.P[i]

        for element in self.elements:
            H = element.H_sum
            ids_matrix = np.zeros((len(element.nodes), len(element.nodes), 2), dtype=int)
            for i in range(len(element.nodes)):
                for j in range(len(element.nodes)):
                    ids_matrix[i][j][0] = element.nodes[i]
                    ids_matrix[j][i][1] = element.nodes[i]
            for i in range(len(ids_matrix)):
                for j in range(len(ids_matrix[i])):
                    H_BC_sum = np.zeros((len(element.nodes), len(element.nodes)))
                    for k in range(len(element.H_BC)):
                        H_BC_sum += element.H_BC[k]
                    self.H_aggregated[ids_matrix[i][j][0]][ids_matrix[i][j][1]] += element.H_sum[i][j] + H_BC_sum[i][j]
                    self.C_aggregated[ids_matrix[i][j][0]][ids_matrix[i][j][1]] += element.C_sum[i][j]


    def min_node(self, nodes, axis):
        minimum = self.nodes[nodes[0]]
        if axis == 'x':
            for i in range(1, len(nodes)):
                if minimum.x > self.nodes[nodes[i]].x:
                    minimum = self.nodes[nodes[i]]
        else:
            for i in range(1, len(nodes)):
                if minimum.y > self.nodes[nodes[i]].y:
                    minimum = self.nodes[nodes[i]]
        return minimum


    def max_node(self, nodes, axis):
        minimum = self.nodes[nodes[0]]
        if axis == 'x':
            for i in range(1, len(nodes)):
                if minimum.x < self.nodes[nodes[i]].x:
                    minimum = self.nodes[nodes[i]]
        else:
            for i in range(1, len(nodes)):
                if minimum.y < self.nodes[nodes[i]].y:
                    minimum = self.nodes[nodes[i]]
        return minimum

        

    def calculate_jacobians(self):
        jacobians = []
        for element in self.elements:
            jacobians.append([])
            for i in range(len(self.gauss_2d_points)):
                pcx = self.deta_derivatives_values[i]
                pcy = self.dksi_derivatives_values[i]
                x = self.getXCoords(element)
                y = self.getYCoords(element)


                result_x = np.sum(pcx*x)
                result_y = np.sum(pcy*y)
                result_xx = np.sum(pcy*x)
                result_yy = np.sum(pcx*y)
                jacobian = np.array([[result_xx, result_x], [result_y, result_yy]])
                jacobians[-1].append(jacobian)
        self.jacobians = np.array(jacobians)
                

    def getXCoords(self, element):
        coords = np.zeros(len(element))
        for i in range(len(element)):
            coords[i] = self.nodes[element[i]].x
        #coords = cartesian_product(coords)
        return coords

    def getYCoords(self, element):
        coords = np.zeros(len(element))
        for i in range(len(element)):
            coords[i] = self.nodes[element[i]].y
        #coords = cartesian_product(coords)
        return coords

    def load_sides(self):
        for i in range(len(self.elements)):
            if self.nodes[self.elements[i].nodes[0]].BC and self.nodes[self.elements[i].nodes[1]].BC:
                points = sort_points(self.gauss_2d_points, True, 1)[:self.integration_points_number]
                nodes = [self.elements[i].nodes[0], self.elements[i].nodes[1]]
                for j in range(len(points)):
                    points[j] = (points[j][0], -1)
                self.elements[i].sides.append(Side(points, nodes, 0, True))
            if self.nodes[self.elements[i].nodes[1]].BC and self.nodes[self.elements[i].nodes[2]].BC:
                points = sort_points(self.gauss_2d_points, False, 0)[:self.integration_points_number]
                nodes = [self.elements[i].nodes[1], self.elements[i].nodes[2]]
                for j in range(len(points)):
                    points[j] = (1, points[j][1])
                self.elements[i].sides.append(Side(points, nodes, 1, False))
            if self.nodes[self.elements[i].nodes[2]].BC and self.nodes[self.elements[i].nodes[3]].BC:
                points = sort_points(self.gauss_2d_points, False, 1)[:self.integration_points_number]
                nodes = [self.elements[i].nodes[2], self.elements[i].nodes[3]]
                for j in range(len(points)):
                    points[j] = (points[j][0], 1)
                self.elements[i].sides.append(Side(points, nodes, 2, True))
            if self.nodes[self.elements[i].nodes[3]].BC and self.nodes[self.elements[i].nodes[0]].BC:
                points = sort_points(self.gauss_2d_points, True, 0)[:self.integration_points_number]
                nodes = [self.elements[i].nodes[3], self.elements[i].nodes[0]]
                for j in range(len(points)):
                    points[j] = (-1, points[j][1])
                self.elements[i].sides.append(Side(points, nodes, 3, False))
                

    shape_funcs = [
                lambda ksi, eta: 0.25*(1-ksi)*(1-eta),
                lambda ksi, eta: 0.25*(1+ksi)*(1-eta),
                lambda ksi, eta: 0.25*(1+ksi)*(1+eta),
                lambda ksi, eta: 0.25*(1-ksi)*(1+eta)
            ]

    shape_funcs_derivatives = [
                [   # dN/deta
                    lambda ksi, eta: (ksi-1)/4,
                    lambda ksi, eta: 0.25*(-ksi-1),
                    lambda ksi, eta: (ksi+1)/4,
                    lambda ksi, eta: (1-ksi)/4
                    ],
                [   # dN/dksi
                    lambda ksi, eta: (eta-1)/4,
                    lambda ksi, eta: (1-eta)/4,
                    lambda ksi, eta: (eta+1)/4,
                    lambda ksi, eta: 0.25*(-eta-1)
                    ]
                ]

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

    gauss_points = []


def cartesian_product(x):
    result = set()
    for first_value in x:
        for second_value in x:
            result.add((first_value, second_value))
    return result

def sides_sort(x, gauss_values, integration_points_number):
        x.sort(key=lambda x: x[1])
        length = len(gauss_values[integration_points_number])
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

def sort_points(points, max, idx):
    cpy = points.copy()
    if max:
        cpy.sort(key=lambda x: x[idx])
    else:
        cpy.sort(key=lambda x: -x[idx])
    return cpy

def pythagoras(x, y):
    return math.sqrt((y.x-x.x)**2+(y.y-x.y)**2)



