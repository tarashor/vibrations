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


def solve(model, mesh, s_matrix, m_matrix):

    s = integrate_matrix(model, mesh, s_matrix)
    m = integrate_matrix(model, mesh, m_matrix)

    fixed_nodes_indicies = mesh.get_fixed_nodes_indicies()

    s = remove_fixed_nodes(s, fixed_nodes_indicies, mesh.nodes_count())
    m = remove_fixed_nodes(m, fixed_nodes_indicies, mesh.nodes_count())


    lam, vec = la.eigh(s, m)

    vec = extend_with_fixed_nodes(vec, fixed_nodes_indicies, mesh.nodes_count())
    
    return lam, vec

def solve_nonlinearity(model, mesh, s_matrix, m_matrix):

    s = integrate_matrix(model, mesh, s_matrix)
    m = integrate_matrix(model, mesh, m_matrix)

    fixed_nodes_indicies = mesh.get_fixed_nodes_indicies()

    s = remove_fixed_nodes(s, fixed_nodes_indicies, mesh.nodes_count())
    m = remove_fixed_nodes(m, fixed_nodes_indicies, mesh.nodes_count())

    lam, vec = la.eigh(s, m)

    res = extend_with_fixed_nodes(vec[:, 0], fixed_nodes_indicies, mesh.nodes_count())
    res_prev = np.zeros(res.shape)

    res = normalize(res)
    lam_nl = lam[0]

    eps = 0.1
    i = 0
    while (np.linalg.norm(res - res_prev) > eps and i < 20):
        res_prev = res
        # print(res_prev.T.shape)
        s_nl = stiffness_nl_matrix(res_prev, model, mesh)
        s_nl = remove_fixed_nodes(s_nl, fixed_nodes_indicies, mesh.nodes_count())
        l_nl, vec_nl = la.eigh(s + s_nl, m)
        print("Freq nl = {}".format(l_nl[0]))
        res = extend_with_fixed_nodes(vec_nl[:, 0], fixed_nodes_indicies, mesh.nodes_count())
        res = normalize(res)
        print("Norm = {}".format(np.linalg.norm(res)))
        lam_nl = l_nl[0]
        i += 1

    return NonlinearResult(lam_nl, res, mesh, model)

def convert_to_results(eigenvalues, eigenvectors, mesh, geometry):
    
    results = []
    for i in range(eigenvalues.size):
        freq = np.sqrt(eigenvalues[i])
        u1 = eigenvectors[:, i][0:mesh.nodes_count()]
        u3 = eigenvectors[:, i][mesh.nodes_count():2 * mesh.nodes_count()]
        u2 = np.zeros((mesh.nodes_count()))
        r = result.Result(freq, u1, u2, u3, mesh, geometry)
        results.append(r)

    return results
    

def integrate_matrix(model, mesh, matrix_func):
    N = 2 * (mesh.nodes_count())
    global_matrix = np.zeros((N, N))
    for element in mesh.elements:
        element_matrix = quadgch5nodes2dim(element_func, element, model.geometry, matrix_func)

        global_matrix += convertToGlobalMatrix(element_matrix, element, N)

    return global_matrix

def integrate_matrix_with_disp(model, mesh, matrix_func, disp):
    N = 2 * (mesh.nodes_count())
    global_matrix = np.zeros((N, N))
    for element in mesh.elements:
        element_matrix = quadgch5nodes2dim(element_func, element, model.geometry, matrix_func, disp)

        global_matrix += convertToGlobalMatrix(element_matrix, element, N)

    return global_matrix

def element_func(ksi, teta, element, geometry, matrix_func):
    x1, x3 = element.to_model_coordinates(ksi, teta)
    x2 = 0
    
    EM = matrix_func(element.material, geometry, x1, x2, x3)
    H = matrices.element_aprox_functions(element, x1, x2, x3)
    J = element.jacobian_element_coordinates()

    e = H.T.dot(EM).dot(H) * J

    return e


def quadgch5nodes2dim(f, element, geometry, matrix_func, disp = None):
    order = 5
    w = [0.23692689, 0.47862867, 0.56888889, 0.47862867, 0.23692689]
    x = [-0.90617985, -0.53846931, 0, 0.53846931, 0.90617985]

    if (disp == None):
        res = w[0] * w[0] * f(x[0], x[0], element, geometry, matrix_func)
    else:
        res = w[0] * w[0] * f(x[0], x[0], element, geometry, matrix_func, disp)

    for i in range(order):
        for j in range(order):
            if (i != 0 or j != 0):
                if (disp == None):
                    res += w[i] * w[j] * f(x[i], x[j], element, geometry, matrix_func)
                else:
                    res += w[i] * w[j] * f(x[i], x[j], element, geometry, matrix_func, disp)

    return res


def map_local_to_global_matrix_index(local_index, element, N):
    global_index = None
    if (local_index % 4 == 0):
        global_index = element.top_left_index
    elif (local_index % 4 == 1):
        global_index = element.top_right_index
    elif (local_index % 4 == 2):
        global_index = element.bottom_right_index
    elif (local_index % 4 == 3):
        global_index = element.bottom_left_index

    if (local_index // 4 == 1):
        global_index += N // 2

    return global_index


def convertToGlobalMatrix(local_matrix, element, N):
    global_matrix = np.zeros((N, N))
    rows, columns = local_matrix.shape
    for i in range(rows):
        for j in range(columns):
            i_global = map_local_to_global_matrix_index(i, element, N)
            j_global = map_local_to_global_matrix_index(j, element, N)
            global_matrix[i_global, j_global] = local_matrix[i, j]

    return global_matrix


def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm
