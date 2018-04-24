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


def solve(model, mesh, t_s_matrix, m_matrix, f_vector, T, time_intervals, u0, v0):

    U = [u0]
    
    eps = 0.001
    
    M = integrate_matrix_with_disp(model, mesh, m_matrix, u0)
    M_inv=  M**-1
    F0 = integrate_vector_with_disp(model, mesh, f_vector, u0)
    
    w0 = -M_inv.dot(F0)
    
    
    for i in range(time_intervals):
        ut = u0
        while True:
            # statement(s)
            
            if not np.linalg.norm(delta_u) < eps:
                break
        L = integrate_vector_with_disp(model, mesh, f_vector, u0)
        delta_u = zeros
        
    K = integrate_matrix_with_disp(model, mesh, t_s_matrix)
    

    fixed_nodes_indicies = mesh.get_fixed_nodes_indicies()

    s = remove_fixed_nodes(s, fixed_nodes_indicies, mesh.nodes_count())
    m = remove_fixed_nodes(m, fixed_nodes_indicies, mesh.nodes_count())


    lam, vec = la.eigh(s, m)

    vec = extend_with_fixed_nodes(vec, fixed_nodes_indicies, mesh.nodes_count())
    
    return lam, vec

def integrate_vector_with_disp(model, mesh, vector_func, disp):
    N = 2 * (mesh.nodes_count())
    global_vector = np.zeros((N, N))
    for element in mesh.elements:
        u_element = matrices.get_u_element(element, disp, mesh.nodes_count())
        element_vector = quadgch5nodes2dim(element_vector_func_disp, element, model.geometry, vector_func, u_element)

        global_vector += convertToGlobalVector(element_vector, element, N)

    return global_vector


def element_vector_func_disp(ksi, teta, element, geometry, vector_func, u_element):
    x1, x3 = element.to_model_coordinates(ksi, teta)
    x2 = 0
    
    grad_u = matrices.get_grad_u(element, u_element, x1, x2, x3)
    
    EF = vector_func(element.material, geometry, x1, x2, x3, grad_u)
    H = matrices.element_aprox_functions(element, x1, x2, x3)
    J = element.jacobian_element_coordinates()

    e = H.T.dot(EF) * J

    return e    

def integrate_matrix_with_disp(model, mesh, matrix_func, disp):
    N = 2 * (mesh.nodes_count())
    global_matrix = np.zeros((N, N))
    for element in mesh.elements:
        u_element = matrices.get_u_element(element, disp, mesh.nodes_count())
        element_matrix = quadgch5nodes2dim(element_func_disp, element, model.geometry, matrix_func, u_element)

        global_matrix += convertToGlobalMatrix(element_matrix, element, N)

    return global_matrix


def element_func_disp(ksi, teta, element, geometry, matrix_func, u_element):
    x1, x3 = element.to_model_coordinates(ksi, teta)
    x2 = 0
    
    grad_u = matrices.get_grad_u(element, u_element, x1, x2, x3)
    
    EM = matrix_func(element.material, geometry, x1, x2, x3, grad_u)
    H = matrices.element_aprox_functions(element, x1, x2, x3)
    J = element.jacobian_element_coordinates()

    e = H.T.dot(EM).dot(H) * J

    return e


def quadgch5nodes2dim(f, element, geometry, matrix_func, disp = None):
    order = 5
    w = [0.23692689, 0.47862867, 0.56888889, 0.47862867, 0.23692689]
    x = [-0.90617985, -0.53846931, 0, 0.53846931, 0.90617985]

    if (disp is None):
        res = w[0] * w[0] * f(x[0], x[0], element, geometry, matrix_func)
    else:
        res = w[0] * w[0] * f(x[0], x[0], element, geometry, matrix_func, disp)

    for i in range(order):
        for j in range(order):
            if (i != 0 or j != 0):
                if (disp is None):
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

def convertToGlobalVector(local_vector, element, N):
    global_vector = np.zeros((N))
    rows, columns = local_vector.shape
    for i in range(rows):
        i_global = map_local_to_global_matrix_index(i, element, N)
        global_vector[i_global] = local_vector[i]

    return global_vector


def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm
