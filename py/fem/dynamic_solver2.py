from . import finiteelements
from . import dmatrices as matrices
from . import result
import numpy as np
# from . import mesh as m
from scipy import linalg as la


def remove_fixed_nodes(matrix, fixed_nodes_indicies, all_nodes_count):
    indicies_to_exclude = i_exclude(fixed_nodes_indicies, all_nodes_count)

    free_nodes1 = [i for i in range(matrix.shape[0]) if i not in indicies_to_exclude]
    free_nodes2 = [i for i in range(matrix.shape[1]) if i not in indicies_to_exclude]
    return matrix[np.ix_(free_nodes1, free_nodes2)]

def remove_fixed_nodes_vector(vector, fixed_nodes_indicies, all_nodes_count):
    indicies_to_exclude = i_exclude(fixed_nodes_indicies, all_nodes_count)

    free_nodes1 = [i for i in range(vector.shape[0]) if i not in indicies_to_exclude]
    return vector[free_nodes1]


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


def solve(model, mesh, t_s_matrix, m_matrix, f_vector, r_vector, T, time_intervals, u0, v0):

    U = [u0]
#    print(u0)
    
    M = integrate_matrix_with_disp(model, mesh, m_matrix, u0)
#    print(np.linalg.norm(u0))
    F0 = integrate_vector_with_disp(model, mesh, f_vector, u0)
    
    R0 = integrate_vector_with_disp(model, mesh, f_vector, u0)
#    Fv0 = integrate_vector_with_disp(model, mesh, f_vector, v0)
#    print(F0-Fv0)
    
    fixed_nodes_indicies = mesh.get_fixed_nodes_indicies()
    
    Mr = remove_fixed_nodes(M, fixed_nodes_indicies, mesh.nodes_count())
    
    F0 = remove_fixed_nodes_vector(F0, fixed_nodes_indicies, mesh.nodes_count())
    R0 = remove_fixed_nodes_vector(R0, fixed_nodes_indicies, mesh.nodes_count())

    
    M_inv= np.linalg.inv(Mr)
#    print(F0.shape)
    
    
    w0 = M_inv.dot(R0-F0)
    
    w0 = extend_with_fixed_nodes(w0, fixed_nodes_indicies, mesh.nodes_count())
#    print(w0)
    
    delta_t = T / time_intervals
    delta_t2 = delta_t * delta_t
    
#    print(M)
    
    M_inv = delta_t2*M_inv
    u_minus_1 = u0-delta_t*v0+delta_t2/2*w0
    
    for i in range(time_intervals):
        print("Iteration = {}".format(i))
        print(u0)
#        print("u0 = {}".format(u0))
        F = integrate_vector_with_disp(model, mesh, f_vector, u0)
        R = integrate_vector_with_disp(model, mesh, r_vector, u0)
        
        R_hat = R-F+2/delta_t2*M.dot(u0)-1/delta_t2*M.dot(u_minus_1)
        
        R_hat = remove_fixed_nodes_vector(R_hat, fixed_nodes_indicies, mesh.nodes_count())
        
        ut = M_inv.dot(R_hat)
        
        ut = extend_with_fixed_nodes(ut, fixed_nodes_indicies, mesh.nodes_count())
        
        u_minus_1 = u0
        u0 = ut
        U.append(ut)
#            
    return U

def integrate_vector_with_disp(model, mesh, vector_func, disp):
    N = 2 * (mesh.nodes_count())
    global_vector = np.zeros((N))
    for element in mesh.elements:
        u_element = matrices.get_u_element(element, disp, mesh.nodes_count())
        element_vector = quadgch5nodes2dim(element_vector_func_disp, element, model.geometry, vector_func, u_element)

        global_vector += convertToGlobalVector(element_vector, element, N)

    return global_vector


def element_vector_func_disp(ksi, teta, element, geometry, vector_func, u_element):
    x1, x3 = element.to_model_coordinates(ksi, teta)
    x2 = 0
    
    grad_u = matrices.get_grad_u(element, geometry, u_element, x1, x2, x3)
    
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
    
    grad_u = matrices.get_grad_u(element, geometry, u_element, x1, x2, x3)
    
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
    rows,  = local_vector.shape
    for i in range(rows):
        i_global = map_local_to_global_matrix_index(i, element, N)
        global_vector[i_global] = local_vector[i]

    return global_vector


def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm
