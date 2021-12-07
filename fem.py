import numpy as np


class Node:
    def __init__(self, x: float = 0.0, y: float = 0.0):
        self.x = x
        self.y = y
        self.BC = False
        self.t = 100.0

    def __repr__(self):
        return f'[{self.x:.3f}, {self.y:.3f}]'


class Element:
    def __init__(self, nodes: np.ndarray = np.zeros(4, dtype=int)):
        self.nodes = nodes
        self.x_derivatives = np.zeros((4,4))
        self.y_derivatives = np.zeros((4,4))
        self.H = np.zeros((4,4,4))
        self.H_sum = np.zeros((4,4))
        self.sides = np.zeros((4,2))
        self.integration_points = np.zeros((4,2,2))
        self.H_BC = np.zeros((4,4,4))
        self.P = np.zeros((4))
        self.C = np.zeros((4,4,4))
        self.C_sum = np.zeros((4,4))


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
    def __init__(self, h: float = 0.0, b: float = 0.0, n_h: float = 0, n_b: float = 0):
        self.c = 700
        self.ro = 7800
        self.H = h
        self.B = b
        self.nH = n_h
        self.nB = n_b
        self.nN = self.nH * self.nB
        self.nE = (self.nH - 1) * (self.nB - 1)
        self.elements = np.zeros((self.nH-1)*(self.nB-1), dtype=object)
        self.nodes = []

        dx = self.B / (self.nB - 1)
        dy = self.H / (self.nH - 1)
        for x in range(self.nB):
            for y in range(self.nH):
                node = Node(float(x * dx), float(y * dy))
                if x == self.nB-1:
                    node.x = self.B
                if y == self.nH-1:
                    node.y = self.H
                if node.x == 0.0 or node.y == 0.0 or node.x == self.B or node.y == self.H:
                    node.BC = True
                self.nodes.append(node)
        self.nodes = np.array(self.nodes)
        '''
        for i in range(self.nH-1):
            for j in range(self.nB-1):
                self.elements[i*(self.nB-1)+j] = Element(np.array([
                    i*self.nH+i, i*self.nH+j+self.nH, i*self.nH+j+self.nH+1, i*self.nH+j+1
                ]))
        '''
        for i in range(self.nB-1):
            for j in range(self.nH-1):
                self.elements[i*(self.nH-1)+j] = Element(np.array([
                    i*(self.nH)+j, (i+1)*(self.nH)+j, (i+1)*(self.nH)+j+1, i*(self.nH)+j+1
                ]))
        self.elements = self.elements.flatten()
        self.load_sides()
        self.H_aggregated = np.zeros((len(self.nodes), len(self.nodes)))
        self.P_aggregated = np.zeros((len(self.nodes)))


    def getXCoords(self, element):
        coords = np.zeros(4)
        for i in range(len(element)):
            coords[i] = self.nodes[element[i]].x
        return coords

    def getYCoords(self, element):
        coords = np.zeros(4)
        for i in range(len(element)):
            coords[i] = self.nodes[element[i]].y
        return coords

    def load_sides(self):
        for i in range(len(self.elements)):
            element = self.elements[i]
            for j in range(len(element)-1):
                if self.nodes[element[j]].BC and self.nodes[element[j+1]].BC:
                    element.sides[j][0] = element[j]
                    element.sides[j][1] = element[j+1]
                    if self.nodes[element[j]].x == self.nodes[element[j+1]]:
                        element.integration_points[j][0][0] = element[j].x
                        element.integration_points[j][0][1] = element[j].y
            if self.nodes[element[-1]].BC and self.nodes[element[0]].BC:
                    element.sides[-1][0] = element[-1]
                    element.sides[-1][1] = element[0]

