import numpy as np


class Node:


    def __init__(self, id, x, y,z):
        self.id = id
        self.x = x
        self.y = y
        self.z = z

        self.index = id - 1  # the index in a vector




class Element:


    def __init__(self, id, node1, node2, node3, node4, node5, node6, node7, node8, node9, node10):
        self.id = id
        self.node1 = node1
        self.node2 = node2
        self.node3 = node3
        self.node4 = node4
        self.node5 = node5
        self.node6 = node6
        self.node7 = node7
        self.node8 = node8
        self.node9 = node9
        self.node10 = node10




class Condition:



    def __init__(self, node, value):
        self.node = node
        self.value = value




class Mesh:

    def __init__(
        self, parameters, nodes, elements, dirichlet_conditions, neumann_conditions
    ):
        self.parameters = np.array(parameters)
        self.nodes = np.array(nodes)
        self.elements = np.array(elements)
        self.dirichlet_conditions = np.array(dirichlet_conditions)
        self.neumann_conditions = np.array(neumann_conditions)


