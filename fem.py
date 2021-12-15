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
    def __init__(self, nodes: np.ndarray, indexes: np.ndarray, horizontal: bool, idx: int):
        self.horizontal = horizontal
        self.nodes = nodes
        self.indexes = indexes
        self.idx = idx


class Element:
    integration_points_number = 2

    def __init__(self, nodes: np.ndarray = np.zeros(4, dtype=int)):
        self.nodes = nodes
        self.x_derivatives = np.zeros((4,4))
        self.y_derivatives = np.zeros((4,4))
        self.H = np.zeros((len(self.nodes), len(self.nodes), len(self.nodes)))
        self.H_sum = np.zeros((len(self.nodes), len(self.nodes)))
        self.integration_points = np.zeros((4,2,2))
        self.H_BC = np.zeros((4, len(self.nodes), len(self.nodes)))
        self.P = np.zeros((len(self.nodes)))
        self.C = np.zeros((len(self.nodes), len(self.nodes), len(self.nodes)))
        self.C_sum = np.zeros((len(self.nodes), len(self.nodes)))
        self.dN_dx = np.zeros((len(self.nodes), len(self.nodes)))
        self.dN_dy = np.zeros((len(self.nodes), len(self.nodes)))
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

    def __init__(self, h: float = 0.0,
            b: float = 0.0,
            n_h: float = 0,
            n_b: float = 0,
            integration_points_number: int = 3):
        Element.integration_points_number = integration_points_number

        gauss_points = [list(cartesian_product(Grid.gauss_values[i].keys()))
                for i in range(len(Grid.gauss_values))]
        for gauss_pts in gauss_points:
            sides_sort(gauss_pts, Grid.gauss_values, integration_points_number-2)



        self.gauss_2d_points = gauss_points[integration_points_number-2]
        self.shape_functions = Grid.shape_funcs[integration_points_number-2]
        self.integration_points_number = integration_points_number

        self.deta_derivatives_values = np.zeros(
                (len(self.gauss_2d_points),
                    len(Grid.shape_funcs[integration_points_number-2]))
                )
        self.dksi_derivatives_values = np.zeros(
                (len(self.gauss_2d_points),
                    len(Grid.shape_funcs[integration_points_number-2]))
                )


        for i in range(len(Grid.shape_funcs_derivatives[integration_points_number-2][0])):
            for j in range(len(self.gauss_2d_points)):
                self.deta_derivatives_values[j][i] = Grid.shape_funcs_derivatives[integration_points_number-2][0][i](*self.gauss_2d_points[j])
                self.dksi_derivatives_values[j][i] = Grid.shape_funcs_derivatives[integration_points_number-2][1][i](*self.gauss_2d_points[j])


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

        dx = self.B / (self.nB+(self.nB-1)*(self.integration_points_number-2)-1)

        cnt = 0
        for x in range((self.nB+(self.nB-1)*(self.integration_points_number-2))):
            number_of_points = 0
            if cnt == 0:
                number_of_points = self.nH+(self.nH-1)*(self.integration_points_number-2)
            else:
                number_of_points = self.nH
            dy = self.H / (number_of_points-1)
            for y in range(number_of_points):
                node = Node(float(x * dx), float(y * dy))
                if x == (self.nB+(self.nB-1)*(self.integration_points_number-2))-1:
                    node.x = self.B
                if y == number_of_points-1:
                    node.y = self.H
                if node.x == 0.0 or node.y == 0.0 or node.x == self.B or node.y == self.H:
                    node.BC = True
                self.nodes.append(node)
            cnt += 1
            if cnt == self.integration_points_number-1:
                cnt = 0

                    
        self.nodes = np.array(self.nodes)

        for i in range(self.nB-1):
            for j in range(self.nH-1):
                idx = []
                # lower side
                for k in range(self.integration_points_number):
                    if k == 0:
                        idx.append(i*((self.integration_points_number-1)*self.nH+(self.nH-1)*(self.integration_points_number-2))+j*(self.integration_points_number-1))
                    elif k == self.integration_points_number-1:
                        idx.append(idx[0]+(self.integration_points_number-1)*self.nH+(self.nH-1)*(self.integration_points_number-2))
                    else:
                        idx.append(idx[0]+self.nH-j*(self.integration_points_number-2)+(self.nH-1)*(self.integration_points_number-2))
                # right side
                for k in range(self.integration_points_number-1):
                    idx.append(idx[-1]+1)
                # upper side
                last_point = idx[-1]
                for k in range(1, self.integration_points_number):
                    if k < self.integration_points_number-1:
                        idx.append(last_point-(k*self.nH+(j+1)*(self.integration_points_number-2)))
                    else:
                        idx.append(idx[-1]-(self.nH+(self.nH-j-2)*(self.integration_points_number-2)))
                # left side
                for k in range(1, self.integration_points_number):
                    if len(idx) < 4*self.integration_points_number - 4:
                        idx.append(idx[-1]-1)
                element = Element(idx)
                self.elements.append(element)


        self.H_aggregated = np.zeros((len(self.nodes), len(self.nodes)))
        self.P_aggregated = np.zeros((len(self.nodes)))
        self.C_aggregated = np.zeros((len(self.nodes), len(self.nodes)))

        self.calculate_jacobians()
        
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


        values_arr = np.zeros((len(self.elements[0].nodes), len(self.elements[0].nodes)))
        for i in range(len(values_arr)):
            for j in range(len(self.shape_functions)):
                values_arr[i][j] = self.shape_functions[j](*self.gauss_2d_points[i]) * Grid.gauss_values[integration_points_number-2][self.gauss_2d_points[i][0]] * Grid.gauss_values[integration_points_number-2][self.gauss_2d_points[i][1]]
        

        # calculate H and C
        for i in range(len(self.elements)):
            for j in range(len(self.jacobians[i])):
                Nx = np.dot(self.elements[i].dN_dx[j,:].reshape(self.elements[i].dN_dx.shape[1], 1),
                        self.elements[i].dN_dx[j,:].reshape(1, -1))
                Ny = np.dot(self.elements[i].dN_dy[j,:].reshape(self.elements[i].dN_dy[j,:].shape[0], 1),
                        self.elements[i].dN_dy[j,:].reshape(1, -1))
                H = self.k * (Nx + Ny) * np.linalg.det(self.jacobians[i][j])
                self.elements[i].H[j] = H
                C = np.zeros((self.integration_points_number**2, self.integration_points_number**2))
                element = self.elements[i]
                element.C_sum += self.c * self.ro * np.dot(values_arr[j].reshape(len(values_arr[j]), 1), values_arr[j].reshape(1, -1)) * np.linalg.det(self.jacobians[i][j])

            for j in range(len(self.jacobians[i])):
                self.elements[i].H_sum += self.elements[i].H[j]

        self.load_sides()


        for element in self.elements:
            for side in element.sides:
                values = np.zeros((self.integration_points_number, len(element.nodes)))
                coords = list(Grid.gauss_values[self.integration_points_number-2].keys())
                for i in range(len(values)):
                    for j in range(len(side.nodes)):
                        if side.idx == 0:
                            values[i][side.indexes[j]] = self.shape_functions[side.indexes[j]](coords[i], -1.0)
                        elif side.idx == 1:
                            values[i][side.indexes[j]] = self.shape_functions[side.indexes[j]](1.0, coords[i])
                        elif side.idx == 2:
                            values[i][side.indexes[j]] = self.shape_functions[side.indexes[j]](coords[i], 1.0)
                        elif side.idx == 3:
                            values[i][side.indexes[j]] = self.shape_functions[side.indexes[j]](-1.0, coords[i])

                axis = 'x'
                if not side.horizontal:
                    axis = 'y'

                minimum_node = self.min_node(side.nodes, axis)
                maximum_node = self.max_node(side.nodes, axis)
                det_J = pythagoras(minimum_node, maximum_node)/self.integration_points_number

                for i in range(len(values)):
                    element.H_BC[side.idx] += self.alpha * (Grid.gauss_values[self.integration_points_number-2][coords[i]]*np.dot(values[i].reshape(len(element.nodes), 1), values[i].reshape(1, -1))) * det_J

        for element in self.elements:
            for side in element.sides:
                values = np.zeros((self.integration_points_number, len(element.nodes)))
                coords = list(Grid.gauss_values[self.integration_points_number-2].keys())
                for i in range(len(values)):
                    for j in range(len(side.nodes)):
                        if side.idx == 0:
                            values[i][side.indexes[j]] = self.shape_functions[side.indexes[j]](coords[i], -1.0)
                        elif side.idx == 1:
                            values[i][side.indexes[j]] = self.shape_functions[side.indexes[j]](1.0, coords[i])
                        elif side.idx == 2:
                            values[i][side.indexes[j]] = self.shape_functions[side.indexes[j]](coords[i], 1.0)
                        elif side.idx == 3:
                            values[i][side.indexes[j]] = self.shape_functions[side.indexes[j]](-1.0, coords[i])

                axis = 'x'
                if not side.horizontal:
                    axis = 'y'

                minimum_node = self.min_node(side.nodes, axis)
                maximum_node = self.max_node(side.nodes, axis)
                det_J = pythagoras(minimum_node, maximum_node)/self.integration_points_number

                for i in range(len(values)):
                    element.P += (self.alpha * (Grid.gauss_values[integration_points_number-2][coords[i]]*values[i].reshape(len(element.nodes), 1)*self.ambient_temperature) * det_J).reshape(len(element.nodes))



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
            for i in range(len(element)):
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
        for element in self.elements:
            nodes = list(element.nodes)

            # vertical sides
            nodes.sort(key=lambda x: self.nodes[x].x)
            # left side
            left_side = nodes[:self.integration_points_number]
            all_true = True
            bc_values = [self.nodes[node].BC for node in left_side]
            for value in bc_values:
                if not value:
                    all_true = False
                    break
            if all_true:
                indexes = []
                for node in left_side:
                    indexes.append(element.nodes.index(node))
                element.sides.append(Side(np.array(left_side), np.array(indexes), False, 3))
            # right side
            right_side = nodes[-self.integration_points_number:]
            all_true = True
            bc_values = [self.nodes[node].BC for node in right_side]
            for value in bc_values:
                if not value:
                    all_true = False
                    break
            if all_true:
                indexes = []
                for node in right_side:
                    indexes.append(element.nodes.index(node))
                element.sides.append(Side(np.array(right_side), np.array(indexes), False, 1))

            # horizontal sides
            nodes.sort(key=lambda x: self.nodes[x].y)
            # upper side
            upper_side = nodes[-self.integration_points_number:]
            all_true = True
            bc_values = [self.nodes[node].BC for node in upper_side]
            for value in bc_values:
                if not value:
                    all_true = False
                    break
            if all_true:
                indexes = []
                for node in upper_side:
                    indexes.append(element.nodes.index(node))
                element.sides.append(Side(np.array(upper_side), np.array(indexes), True, 2))
            # lower side
            lower_side = nodes[:self.integration_points_number]
            all_true = True
            bc_values = [self.nodes[node].BC for node in lower_side]
            for value in bc_values:
                if not value:
                    all_true = False
                    break
            if all_true:
                indexes = []
                for node in lower_side:
                    indexes.append(element.nodes.index(node))
                element.sides.append(Side(np.array(lower_side), np.array(indexes), True, 0))


    shape_funcs = [
            [
                lambda ksi, eta: 0.25*(1-ksi)*(1-eta),
                lambda ksi, eta: 0.25*(1+ksi)*(1-eta),
                lambda ksi, eta: 0.25*(1+ksi)*(1+eta),
                lambda ksi, eta: 0.25*(1-ksi)*(1+eta)
                ],
            [
                lambda ksi, eta: -0.25*(1-ksi)*(1-eta)*(ksi+eta+1),
                lambda ksi, eta: 0.5*(1-ksi**2)*(1-eta),
                lambda ksi, eta: 0.25*(1+ksi)*(1-eta)*(ksi-eta-1),
                lambda ksi, eta: 0.5*(1-eta**2)*(1+ksi),
                lambda ksi, eta: 0.25*(1+ksi)*(1+eta)*(ksi+eta-1),
                lambda ksi, eta: 0.5*(1-ksi**2)*(1+eta),
                lambda ksi, eta: -0.25*(1-ksi)*(1+eta)*(ksi-eta+1),
                lambda ksi, eta: 0.5*(1-eta**2)*(1-ksi)
                ]
            ]

    shape_funcs_derivatives = [
            [
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
                ],
            [
                [   # dN/deta
                    lambda ksi, eta: -0.25*(ksi-1)*(2*eta+ksi),
                    lambda ksi, eta: 0.5*ksi**2-0.5,
                    lambda ksi, eta: 0.25*(ksi+1)*(2*eta-ksi),
                    lambda ksi, eta: -eta*(ksi+1),
                    lambda ksi, eta: 0.25*(ksi+1)*(2*eta+ksi),
                    lambda ksi, eta: 0.5-0.5*ksi**2,
                    lambda ksi, eta: -0.25*(ksi-1)*(2*eta-ksi),
                    lambda ksi, eta: eta*(ksi-1)
                    ],
                [   # dN/dksi
                    lambda ksi, eta: -0.25*(eta-1)*(eta+2*ksi),
                    lambda ksi, eta: (eta-1)*ksi,
                    lambda ksi, eta: 0.25*(eta-1)*(eta-2*ksi),
                    lambda ksi, eta: 0.5*(1-eta**2),
                    lambda ksi, eta: 0.25*(eta+1)*(eta+2*ksi),
                    lambda ksi, eta: -(eta+1)*ksi,
                    lambda ksi, eta: -0.25*(eta+1)*(eta-2*ksi),
                    lambda ksi, eta: 0.5*eta**2-0.5
                    ]
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

def pythagoras(x, y):
    return math.sqrt((y.x-x.x)**2+(y.y-x.y)**2)



