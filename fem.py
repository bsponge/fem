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

    def __repr__(self):
        return self.nodes.__repr__()


class Grid:
    def __init__(self, h: float = 0.0, b: float = 0.0, n_h: float = 0, n_b: float = 0):
        self.H = h
        self.B = b
        self.nH = n_h
        self.nB = n_b
        self.nN = self.nH * self.nB
        self.nE = (self.nH - 1) * (self.nB - 1)
        self.elements = np.empty(shape=(self.nH, self.nB), dtype=object)
        self.nodes = np.zeros(self.nN, dtype=object)
        for x in range(self.nH):
            for y in range(self.nB):
                element = Element()
                self.elements[x][y] = element
        dx = self.H / (self.nH - 1)
        dy = self.B / (self.nB - 1)
        for x in range(self.nH):
            for y in range(self.nB):
                node = Node(float(x * dx), float(y * dy))
                self.nodes[x * self.nB + y] = node
        for x in range(self.elements.shape[0]):
            for y in range(self.elements.shape[1]):
                self.elements[x][y].nodes = np.array([y + x * self.nH,
                                                      y + self.nH + x * self.nH,
                                                      y + self.nH + x * self.nH + 1,
                                                      y + x * self.nH + 1])
        self.elements = self.elements.flatten()



