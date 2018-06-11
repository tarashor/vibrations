from . import matrices1D as matrices
from ...model import Model as mod
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
    fixed_indicies3 = [3 * x + 2 for x in fixed_nodes_indicies]
    fixed_indicies1 = [3 * x for x in fixed_nodes_indicies]
    fixed_indicies2 = [3 * x + 1 for x in fixed_nodes_indicies]
    
#    fixed_indicies1 = []
#    fixed_indicies2 = []
    return sorted(fixed_indicies1+fixed_indicies2+fixed_indicies3)


def solve(model, mesh, s_matrix, m_matrix):

    s = integrate_matrix(model, mesh, s_matrix)
    
    m = integrate_matrix(model, mesh, m_matrix)
    

    fixed_nodes_indicies = mesh.get_fixed_nodes_indicies()

    s = remove_fixed_nodes(s, fixed_nodes_indicies, mesh.nodes_count())
    m = remove_fixed_nodes(m, fixed_nodes_indicies, mesh.nodes_count())
    
#    print(m)
    
#    if (model.boundary_conditions == mod.FIXED_BOTTOM_LEFT_RIGHT_POINTS):
#        cols = s.shape[1]
#        
#        s_h = np.zeros((2, cols))
#        s_h[0,0] = 1
#        s_h[0,2] = model.fixed_x3
#        s_h[1,cols-3] = 1
#        s_h[1,cols-2] = model.fixed_x3
#        
#        s = np.concatenate((s, s_h), axis=0)
#        
#        s_h = np.insert(s_h, s_h.shape[1], 0, axis=1)
#        s_h = np.insert(s_h, s_h.shape[1], 0, axis=1)
#        
#        s = np.concatenate((s, s_h.T), axis=1)
#        
#        m = np.insert(m, m.shape[0], 0, axis=0)
#        m = np.insert(m, m.shape[0], 0, axis=0)
#        m = np.insert(m, m.shape[1], 0, axis=1)
#        m = np.insert(m, m.shape[1], 0, axis=1)
        
#    print(m)
#    evm = la.eigvalsh(m)
#    print(evm)

#    lam, vec = la.eig(s, m, left=False, right=True)
    lam, vec = la.eigh(s, m)
    
#    print(lam)

    vec = extend_with_fixed_nodes(vec, fixed_nodes_indicies, mesh.nodes_count())
    
    return lam, vec

    

def integrate_matrix(model, mesh, matrix_func):
    N = 3 * (mesh.nodes_count())
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
    if (local_index // 3 == 0):
        global_index = 3*element.start_index+(local_index % 3)
    else:
        global_index = 3*element.end_index+(local_index % 3)

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

