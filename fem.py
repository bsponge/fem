import numpy as np


class Node:
    def __init__(self, x: float = 0.0, y: float = 0.0):
        self.x = x
        self.y = y

    def __repr__(self):
        return f'[{self.x:.3f}, {self.y:.3f}]'


class Element:
    def __init__(self, nodes: np.ndarray = np.zeros(4, dtype=int)):
        self.nodes = nodes
        self.x_derivatives = np.zeros((4,4))
        self.y_derivatives = np.zeros((4,4))
        self.H = np.zeros((4,4,4))
        self.H_sum = np.zeros((4,4))

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
        self.H = h
        self.B = b
        self.nH = n_h
        self.nB = n_b
        self.nN = self.nH * self.nB
        self.nE = (self.nH - 1) * (self.nB - 1)
        self.elements = np.zeros((self.nH-1)*(self.nB-1), dtype=object)
        self.nodes = []#np.zeros(self.nN, dtype=object)
        dx = self.B / (self.nB - 1)
        dy = self.H / (self.nH - 1)
        for x in range(self.nB):
            for y in range(self.nH):
                node = Node(float(x * dx), float(y * dy))
                self.nodes.append(node)
        self.nodes = np.array(self.nodes)
        '''
        for x in range(self.elements.shape[0]):
            for y in range(self.elements.shape[1]):
                self.elements[x][y].nodes = np.array([y + x * self.nH,
                                                      y + self.nH + x * self.nH,
                                                      y + self.nH + x * self.nH + 1,
                                                      y + x * self.nH + 1])
                                                      '''
        for i in range(self.nB-1):
            for j in range(self.nH-1):
                self.elements[i*self.nB+j] = Element(np.array([i*self.nH+j, i*self.nH+j+self.nH, i*self.nH+j+self.nH+1, i*self.nH+j+1]))

        self.elements = self.elements.flatten()


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



