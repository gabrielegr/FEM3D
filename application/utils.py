import numpy as np
from functools import reduce
from .classes import Node, Element, Mesh, Condition


def find_node_in_nodes(node_id, nodes):

    return next((node for node in nodes if node.id == node_id))


def read_mesh(filename):
        with open(f"data/{filename}.dat") as f_in:

            lines = list(
            map(
                lambda line: line.split(" "),
                filter(None, (line.rstrip() for line in f_in)),
            )
            )

        lines = [list(filter(lambda x: x != "", line)) for line in lines]

        current_line = 0

        E1, F1, F2, F3 = map(float, lines[current_line])
        current_line += 1

        cant_nodes, cant_elements, cant_dirichletX, cant_dirichletY, cant_dirichletZ, cant_neumann = map(int, lines[current_line])


        current_line += 2

        nodes = []
        for line in lines[current_line : current_line + cant_nodes]:

            node_id, node_x, node_y, node_z = line


            node = Node(int(node_id), float(node_x), float(node_y), float(node_z))


            nodes.append(node)


            current_line += 1


        current_line += 2
        elements = []
        for line in lines[current_line : current_line + cant_elements]:

            element_id, element_node1_id, element_node2_id, element_node3_id, element_node4_id, element_node5_id, element_node6_id,  element_node7_id,  element_node8_id,  element_node9_id,  element_node10_id = map(
                int, line
            )


            element_node1 = find_node_in_nodes(element_node1_id, nodes)
            element_node2 = find_node_in_nodes(element_node2_id, nodes)
            element_node3 = find_node_in_nodes(element_node3_id, nodes)
            element_node4 = find_node_in_nodes(element_node4_id, nodes)
            element_node5 = find_node_in_nodes(element_node5_id, nodes)
            element_node6 = find_node_in_nodes(element_node6_id, nodes)
            element_node7 = find_node_in_nodes(element_node7_id, nodes)
            element_node8 = find_node_in_nodes(element_node8_id, nodes)
            element_node9 = find_node_in_nodes(element_node9_id, nodes)
            element_node10 = find_node_in_nodes(element_node10_id, nodes)

            element = Element(element_id, element_node1, element_node2, element_node3, element_node4, element_node5, element_node6, element_node7, element_node8, element_node9, element_node10)


            elements.append(element)


            current_line += 1


        current_line += 2

        dirichletX_conditions = []
        for line in lines[current_line : current_line + cant_dirichletX]:

            condition_node_id, condition_value = line

            condition_node = find_node_in_nodes(int(condition_node_id), nodes)

            condition = Condition(condition_node, float(condition_value))

            dirichletX_conditions.append(condition)


            current_line += 1


        current_line += 2
        dirichletY_conditions = []
        for line in lines[current_line : current_line + cant_dirichletY]:

            condition_node_id, condition_value = line


            condition_node = find_node_in_nodes(int(condition_node_id), nodes)


            condition = Condition(condition_node, float(condition_value))

            dirichletY_conditions.append(condition)


            current_line += 1

        current_line += 2
        dirichletZ_conditions = []
        for line in lines[current_line : current_line + cant_dirichletZ]:

            condition_node_id, condition_value = line

            condition_node = find_node_in_nodes(int(condition_node_id), nodes)

            condition = Condition(condition_node, float(condition_value))

            dirichletZ_conditions.append(condition)


            current_line += 1

        dirichlet_conditions = []    
        for i in range(114):
            dirichlet_conditions.append([dirichletX_conditions[i],dirichletY_conditions[i],dirichletZ_conditions[i]])

        current_line += 2
        neumann_conditions = []
        for line in lines[current_line : current_line + cant_neumann]:

            condition_node_id, condition_value = line

            condition_node = find_node_in_nodes(int(condition_node_id), nodes)


            condition = Condition(condition_node, float(condition_value))

            neumann_conditions.append(condition)

        return Mesh([E1, [F1, F2, F3]], nodes, elements, dirichlet_conditions, neumann_conditions)


def assembly_local_K(element, EI):
    

    node1 = element.node1
    node2 = element.node2
    node3 = element.node3
    node4 = element.node4
    node5 = element.node5
    node6 = element.node6
    node7 = element.node7
    node8 = element.node8
    node9 = element.node9
    node10 = element.node10

    c1 = 1 / pow(0.0000001 - node1.x, 2)
    c2 = (4 * node1.x + 4 * node2.x - 8 * node8.x)/(0.0000001 - node1.x)

    if c1 == 0:
        c1 = 0.0000001
    if c2 == 0:
        c2 = 0.0000001

    Ja = np.linalg.det(
    np.array([
        [(node2.x - node1.x), (node3.x - node1.x), (node4.x - node1.x)],
        [(node2.y - node1.y), (node3.y - node1.y), (node4.y - node1.y)],
        [(node2.z - node1.z), (node3.z - node1.z), (node4.z - node1.z)],
            ])
    )


    A = -(pow(4 * c1 - c2, 4) / (192 * pow(c2, 2))) - (pow(4 * c1 - c2, 3) / (24 * c2)) - (pow(4 * c1 - c2, 5) / (3840 * pow(c2, 3))) + (pow(4 * c1 + 3 * c2, 5) / (3840 * pow(c2, 3)))

    B = -(pow(4 * c1 + c2, 4) / (192 * pow(c2, 2))) + (pow(4 * c1 + c2, 3) / (24 * c2)) + (pow(4 * c1 + c2, 5) / (3840 * pow(c2, 3))) - (pow(4 * c1 - 3 * c2, 5) / (3840 * pow(c2, 3)))

    C = 4/15 * pow(c2,2) 

    D = (pow(4 * c2 - c1, 4)/(192 * pow(c2, 2))) - (pow(4 * c2 - c1, 5)/(3840 * pow(c2, 3))) + (pow(4 * c2 + 8 * c1, 5) / (7680 * pow(c2, 3))) - (pow(4 * c2 - 8 * c1, 5) * (7 / (7680 * pow(c2, 3)))) + (pow(-8 * c1, 5) / (768 * pow(c2, 3))) - (pow(4 * c2 - 8 * c1, 4) * (c1 / (96 * pow(c2, 3)))) + ((pow(-8 * c1, 4)) * ((2 * c1 - 1) / (192 * pow(c2, 3))))

    E = (8/3 * pow(c1, 2)) + (1/30 * pow(c2, 2)) 

    F = (2/3 * c1 * c2) - (1/30 * pow(c2, 2))

    G = -(16/3 * pow(c1, 2)) - (4/3 * c1 * c2) - (2/15 * pow(c2, 2))

    H = (2/3 * c1 * c2) + (1/30 * pow(c2, 2))

    I = -(16/3 * pow(c1, 2)) - (2/3 * pow(c2, 2))

    J = (2/15 * pow(c2, 2))

    K = (-4/3 * c1 * c2)

    U = np.array([
    [A, E, 0, 0, -1*F, 0, -1*F, G, F, F],
    [E, B, 0, 0, -1*H, 0, -1*H, I, H, H],
    [0, 0, 0, 0, 0,  0,  0, 0, 0, 0],
    [0, 0, 0, 0, 0,  0,  0, 0, 0, 0],
    [-1*F, -1*H, 0, 0, C, 0, J, -1*K, -1*C, -1*J],
    [0, 0, 0, 0, 0,  0,  0, 0, 0, 0],
    [-1*F, -1*H, 0, 0, J, 0, C, -1*K, -1*J, -1*C],
    [G, I, 0, 0, -1*K,  0,  -1*K, D, K, K],
    [F, H, 0, 0, -1*C,  0,  -1*J, K, C, J],
    [F, H, 0, 0, -1*J,  0,  -1*C, K, J, C],
    ])
    uz= np.zeros((10, 10))
    Um =  np.array([
    [U, uz, uz],
    [uz, U, uz],
    [uz, uz, U],
    ])

    Kl = EI * Ja * Um

    return Kl


def assembly_local_b(element, F):

    node1 = element.node1
    node2 = element.node2
    node3 = element.node3
    node4 = element.node4
    node5 = element.node5
    node6 = element.node6
    node7 = element.node7
    node8 = element.node8
    node9 = element.node9
    node10 = element.node10


    Ja = np.linalg.det(
    np.array(
    [
        [(node2.x - node1.x), (node3.x - node1.x), (node4.x - node1.x)],
        [(node2.y - node1.y), (node3.y - node1.y), (node4.y - node1.y)],
        [(node2.z - node1.z), (node3.z - node1.z), (node4.z - node1.z)],
    ]
       )
    )


    T = np.array([[59], [-1], [-1], [-1], [4], [4], [4], [4], [4], [4]])

    tz = np.zeros((10, 10))
    Tm =  np.array([
    [T, tz, tz],
    [tz, T, tz],
    [tz, tz, T],
    ]
    )

    Bl = (Ja/120) * Tm * F

    return Bl


def assembly(nodes,elements, localKs, localbs):

    K = np.zeros((178, 178))
    b = np.zeros(178)

    for index in range(178):
        b[index] = 0
        for count in range(178):
            K[index][count] = 0

    for index, element in enumerate(elements):
        # get nodes from element
        node1 = element.node1
        node2 = element.node2
        node3 = element.node3
        node4 = element.node4
        node5 = element.node5
        node6 = element.node6
        node7 = element.node7
        node8 = element.node8
        node9 = element.node9
        node10 = element.node10

        localK = localKs[index]

        K[node1.index][node1.index] = np.reduce(localK[0][0], K[node1.index][node1.index])
        K[node1.index][node2.index] = np.reduce(localK[0][1], K[node1.index][node2.index])
        K[node1.index][node3.index] = np.reduce(localK[0][2], K[node1.index][node3.index])

        K[node2.index][node1.index] = np.reduce(localK[1][0], K[node1.index][node1.index])
        K[node2.index][node2.index] = np.reduce(localK[1][1], K[node1.index][node2.index])
        K[node2.index][node3.index] = np.reduce(localK[1][2], K[node1.index][node3.index])

        K[node3.index][node1.index] = np.reduce(localK[2][0], K[node1.index][node1.index])
        K[node3.index][node2.index] = np.reduce(localK[2][1], K[node1.index][node2.index])
        K[node3.index][node3.index] = np.reduce(localK[2][2], K[node1.index][node3.index])

        K[node4.index][node1.index] = np.reduce(localK[3][0], K[node1.index][node1.index])
        K[node4.index][node2.index] = np.reduce(localK[3][1], K[node1.index][node2.index])
        K[node4.index][node3.index] = np.reduce(localK[3][2], K[node1.index][node3.index])

        K[node5.index][node1.index] = np.reduce(localK[4][0], K[node1.index][node1.index])
        K[node5.index][node2.index] = np.reduce(localK[4][1], K[node1.index][node2.index])
        K[node5.index][node3.index] = np.reduce(localK[4][2], K[node1.index][node3.index])

        K[node6.index][node1.index] = np.reduce(localK[5][0], K[node1.index][node1.index])
        K[node6.index][node2.index] = np.reduce(localK[5][1], K[node1.index][node2.index])
        K[node6.index][node3.index] = np.reduce(localK[5][2], K[node1.index][node3.index])

        K[node7.index][node1.index] = np.reduce(localK[6][0], K[node1.index][node1.index])
        K[node7.index][node2.index] = np.reduce(localK[6][1], K[node1.index][node2.index])
        K[node7.index][node3.index] = np.reduce(localK[6][2], K[node1.index][node3.index])

        K[node8.index][node1.index] = np.reduce(localK[7][0], K[node1.index][node1.index])
        K[node8.index][node2.index] = np.reduce(localK[7][1], K[node1.index][node2.index])
        K[node8.index][node3.index] = np.reduce(localK[7][2], K[node1.index][node3.index])

        K[node9.index][node1.index] = np.reduce(localK[8][0], K[node1.index][node1.index])
        K[node9.index][node2.index] = np.reduce(localK[8][1], K[node1.index][node2.index])
        K[node9.index][node3.index] = np.reduce(localK[8][2], K[node1.index][node3.index])

        K[node10.index][node1.index] = np.reduce(localK[9][0], K[node1.index][node1.index])
        K[node10.index][node2.index] = np.reduce(localK[9][1], K[node1.index][node2.index])
        K[node10.index][node3.index] = np.reduce(localK[9][2], K[node1.index][node3.index])



        local_b = localbs[index]
        b[node1.index] += local_b[0]
        b[node2.index] += local_b[1]
        b[node3.index] += local_b[2]
        b[node4.index] += local_b[3]

    return K, b


def apply_conditions(neumann_conditions, dirichlet_conditions, K, b):

    temp_K = K
    temp_b = b


    for neumann_condition in neumann_conditions:
        neumann_node = neumann_condition.node
        neumann_value = neumann_condition.value
        temp_b[neumann_node.index] += neumann_value


    sorted_conditions = sorted(
        dirichlet_conditions, key=lambda condition: condition.node.id
    )


    for i, dirichlet_condition in enumerate(sorted_conditions):
        dirichlet_node = dirichlet_condition.node
        dirichlet_value = dirichlet_condition.value

        temp_K = np.delete(temp_K, dirichlet_node.index - i, axis=0)

        temp_b = np.delete(temp_b, dirichlet_node.index - i)

        for index, row in enumerate(temp_K):
            cell = row[dirichlet_node.index - i]
            temp_b[index] += -1 * dirichlet_value * cell

        # delete from K, the object (index) from axis 1 (column)
        temp_K = np.delete(temp_K, dirichlet_node.index - i, axis=1)

    return temp_K, temp_b


def calculate_fem(K, b, conditions):

    K_inv = np.linalg.inv(K)
    X = K_inv.dot(b)

    for condition in conditions:
        X = np.insert(X, condition.node.index, condition.value)

    return list(map(lambda x: round(x, 2), X))

def post_processing_input(filename, X):

    with open(f"data/{filename}.post.res", "w") as f_out:
        f_out.write("Results\n")
        f_out.write('Result "Movement" "Load Case 1" 1 Scalar OnNodes\n')
        f_out.write('ComponentNames "X"\n')
        f_out.write("Values\n")

        for i, value in enumerate(X):
            f_out.write("{0}\t{1}\n".format(i + 1, value))
        f_out.write("End Values")