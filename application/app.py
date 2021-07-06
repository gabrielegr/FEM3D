import sys
import numpy as np

from .utils import (
    read_mesh,
    assembly_local_K,
    assembly_local_b,
    assembly,
    apply_conditions,
    calculate_fem,
    post_processing_input
)


def run():

    args = sys.argv[1:]

    data file_name = 'viga'


    mesh = read_mesh(data_filename)


    EI, F = mesh.parameters
    elements = mesh.elements
    nodes = mesh.nodes


    local_K = np.array([assembly_local_K(element, EI) for element in elements])
    local_b= np.array([assembly_local_b(element, F) for element in elements])

    print(local_K[0])


    K, b = assembly(nodes, elements, local_K, local_b)


    K, b = apply_conditions(mesh.neumann_conditions, mesh.dirichlet_conditions, K, b)

    X = calculate_fem(K, b, mesh.dirichlet_conditions)

    print("X: {0}".format(X))
    post_processing_input(data_filename, X)
