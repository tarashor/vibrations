from . import finiteelements
from . import matrices
from . import result
import numpy as np
# from . import mesh as m
from scipy import linalg as la

def remove_fixed_nodes(matrix, fixed_nodes_indicies, all_nodes_count):
    indicies_to_exclude = i_exclude(fixed_nodes_indicies, all_nodes_count)

    free_nodes1 = [i for i in range(matrix.shape[0]) if i not in indicies_to_exclude]
    free_nodes2 = [i for i in range(matrix.shape[1]) if i not in indicies_to_exclude]
    return matrix[np.ix_(free_nodes1, free_nodes2)]


def extend_with_fixed_nodes(eig_vectors, fixed_nodes_indicies, all_nodes_count):
    indicies_to_exclude = i_exclude(fixed_nodes_indicies, all_nodes_count)
    res = eig_vectors
    for i in indicies_to_exclude:
        res = np.insert(res, i, 0, axis=0)

    return res


def i_exclude(fixed_nodes_indicies, nodes_count):
    fixed_u1_indicies = fixed_nodes_indicies
    fixed_u2_indicies = [nodes_count + x for x in fixed_nodes_indicies]
    return sorted(fixed_u1_indicies + fixed_u2_indicies)


def solve(model, mesh):

    s = stiffness_matrix(model, mesh)
    m = mass_matrix(model, mesh)

    fixed_nodes_indicies = mesh.get_fixed_nodes_indicies()

    s = remove_fixed_nodes(s, fixed_nodes_indicies, mesh.nodes_count())
    m = remove_fixed_nodes(m, fixed_nodes_indicies, mesh.nodes_count())
    
#    print(s)
#    print(m)

    lam, vec = la.eigh(s, m)

    vec = extend_with_fixed_nodes(vec, fixed_nodes_indicies, mesh.nodes_count())
    
    results = []
    for i in range(4):
        freq = np.sqrt(lam[i])
        u1 = vec[:, i][0:mesh.nodes_count()]
        u2 = vec[:, i][mesh.nodes_count():2 * mesh.nodes_count()]
        u3 = np.zeros((mesh.nodes_count()))
        r = result.Result(freq, u1, u2, u3, mesh, model.geometry)
        results.append(r)
    

    return results


def stiffness_matrix(model, mesh):
    N = 2 * (mesh.nodes_count())
    K = np.zeros((N, N))
    for element in mesh.elements:
        K_element = quadgch5nodes2dim(k_element_func, element, model.geometry)

        K += convertToGlobalMatrix(K_element, element, N)

    return K



def k_element_func(ksi, teta, element, geometry):
    x1, x2 = element.to_model_coordinates(ksi, teta)
    x3 = 0
#    print("ksi = {}, teta = {}".format(ksi, teta))
#    print("alpha1 = {}, alpha2 = {}".format(x1, x2))
    C = element.material.tensor_C(geometry, x1, x2, x3)
    E = matrices.grad_to_strain()
    B = matrices.deriv_to_grad(geometry, x1, x2, x3)
#    print("B={}".format(B))
    N = matrices.element_aprox_functions(element, x1, x2, x3)
#    print("N={}".format(N))
    J = element.jacobian_element_coordinates()
    
    g = geometry.metric_tensor(x1, x2, x3)
    
    k=N.T.dot(B.T).dot(E.T).dot(C).dot(E).dot(B).dot(N) * J * np.linalg.det(g)
#    print("k={}".format(k))

    return k


def mass_matrix(model, mesh):
    N = 2 * (mesh.nodes_count())
    M = np.zeros((N, N))
    for element in mesh.elements:
        M_element = quadgch5nodes2dim(m_element_func, element, model.geometry)
        M += convertToGlobalMatrix(M_element, element, N)

    return M

def m_element_func(ksi, teta, element, geometry):
    x1, x2 = element.to_model_coordinates(ksi, teta)
    x3 = 0
    g = geometry.metric_tensor(x1, x2, x3)
    B_s = matrices.deriv_to_vect()
    N = matrices.element_aprox_functions(element, x1, x2, x3)
    J = element.jacobian_element_coordinates()
    return element.material.rho * N.T.dot(B_s.T.dot(g.dot(B_s.dot(N)))) * J * np.linalg.det(g)


def quadgch5nodes2dim(f, element, geometry):
    order=5
    w=[0.23692689, 0.47862867, 0.56888889, 0.47862867, 0.23692689]
    x=[-0.90617985, -0.53846931, 0, 0.53846931, 0.90617985]

    res=w[0] * w[0] * f(x[0], x[0], element, geometry)

    for i in range(order):
        for j in range(order):
            if (i != 0 or j != 0):
                res += w[i] * w[j] * f(x[i], x[j], element, geometry)

    return res


def map_local_to_global_matrix_index(local_index, element, N):
    global_index=None
    if (local_index % 4 == 0):
        global_index=element.top_left_index
    elif (local_index % 4 == 1):
        global_index=element.top_right_index
    elif (local_index % 4 == 2):
        global_index=element.bottom_right_index
    elif (local_index % 4 == 3):
        global_index=element.bottom_left_index

    if (local_index // 4 == 1):
        global_index += N // 2

    return global_index


def convertToGlobalMatrix(local_matrix, element, N):
    global_matrix=np.zeros((N, N))
    rows, columns=local_matrix.shape
    for i in range(rows):
        for j in range(columns):
            i_global=map_local_to_global_matrix_index(i, element, N)
            j_global=map_local_to_global_matrix_index(j, element, N)
            global_matrix[i_global, j_global]=local_matrix[i, j]

    return global_matrix


def normalize(v):
    norm=np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm
