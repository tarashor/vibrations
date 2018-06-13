from . import matrices1D as matrices
from ...model import Model as mod
import numpy as np
# from . import mesh as m
from scipy import linalg as la


def remove_fixed_nodes(matrix, fixed_nodes_indicies, all_nodes_count, bc):
    indicies_to_exclude = i_exclude(fixed_nodes_indicies, all_nodes_count, bc)

    free_nodes1 = [i for i in range(matrix.shape[0]) if i not in indicies_to_exclude]
    free_nodes2 = [i for i in range(matrix.shape[1]) if i not in indicies_to_exclude]
    return matrix[np.ix_(free_nodes1, free_nodes2)]


def extend_with_fixed_nodes(eig_vectors, fixed_nodes_indicies, all_nodes_count, bc):
    indicies_to_exclude = i_exclude(fixed_nodes_indicies, all_nodes_count, bc)
    res = eig_vectors
    for i in indicies_to_exclude:
        res = np.insert(res, i, 0, axis=0)

    return res


def i_exclude(fixed_nodes_indicies, nodes_count, bound_cond):
    fixed_indicies1= [6 * x for x in fixed_nodes_indicies]
    fixed_indicies4 = [6 * x + 3 for x in fixed_nodes_indicies]
    
    fixed_indicies2 = []
    fixed_indicies3 = []
    fixed_indicies5 = []
    fixed_indicies6 = []
    
    if (bound_cond == mod.FIXED_LEFT_RIGHT_EDGE):
        fixed_indicies2 = [6 * x + 1 for x in fixed_nodes_indicies]
        fixed_indicies3 = [6 * x + 2 for x in fixed_nodes_indicies]    
        fixed_indicies5 = [6 * x + 4 for x in fixed_nodes_indicies]    
        fixed_indicies6 = [6 * x + 5 for x in fixed_nodes_indicies]    
    
    return sorted(fixed_indicies1+fixed_indicies2+
                  fixed_indicies3+fixed_indicies4+
                  fixed_indicies5+fixed_indicies6)


def solve(model, mesh, s_matrix, m_matrix):

    s = integrate_matrix(model, mesh, s_matrix)
    
    m = integrate_matrix(model, mesh, m_matrix)
    

    fixed_nodes_indicies = mesh.get_fixed_nodes_indicies()

    s = remove_fixed_nodes(s, fixed_nodes_indicies, mesh.nodes_count(), model.boundary_conditions)
    m = remove_fixed_nodes(m, fixed_nodes_indicies, mesh.nodes_count(), model.boundary_conditions)

    lam, vec = la.eigh(s, m)

    vec = extend_with_fixed_nodes(vec, fixed_nodes_indicies, mesh.nodes_count(), model.boundary_conditions)
    
    return lam, vec

    

def integrate_matrix(model, mesh, matrix_func):
    N = 6 * (mesh.nodes_count())
    global_matrix = np.zeros((N, N))
    for element in mesh.elements:
        element_matrix = quadgch5nodes1dim(element_func, element, model.geometry, matrix_func)

        global_matrix += convertToGlobalMatrix(element_matrix, element, N)

    return global_matrix

def element_func(ksi, element, geometry, matrix_func):
    x1 = element.to_model_coordinates(ksi)
    x3 = 0
    x2 = 0
    
    EM = matrix_func(element.material, geometry, x1, element.thickness)
    H = matrices.element_aprox_functions(element, x1, x2, x3)
    J = element.jacobian_element_coordinates()

    e = H.T.dot(EM).dot(H) * J

    return e


def quadgch5nodes1dim(f, element, geometry, matrix_func, disp = None):
    order = 5
    w = [0.23692689, 0.47862867, 0.56888889, 0.47862867, 0.23692689]
    x = [-0.90617985, -0.53846931, 0, 0.53846931, 0.90617985]

    if (disp is None):
        res = w[0] * f(x[0], element, geometry, matrix_func)
    else:
        res = w[0] * f(x[0], element, geometry, matrix_func, disp)

    for i in range(order):
        if (i != 0):
            if (disp is None):
                res += w[i] * f(x[i], element, geometry, matrix_func)
            else:
                res += w[i] * f(x[i], element, geometry, matrix_func, disp)

    return res


def map_local_to_global_matrix_index(local_index, element, N):
    global_index = None
    if (local_index // 6 == 0):
        global_index = 6*element.start_index+(local_index % 6)
    else:
        global_index = 6*element.end_index+(local_index % 6)

    return global_index


def convertToGlobalMatrix(local_matrix, element, N):
    global_matrix = np.zeros((N, N))
    rows, columns = local_matrix.shape
    for i in range(rows):
        for j in range(columns):
#            print(element)
            i_global = map_local_to_global_matrix_index(i, element, N)
            j_global = map_local_to_global_matrix_index(j, element, N)
            
#            print('{} - {}'.format(i, i_global))
            
#            print('{} - {}'.format(j, j_global))
            global_matrix[i_global, j_global] = local_matrix[i, j]

    return global_matrix

