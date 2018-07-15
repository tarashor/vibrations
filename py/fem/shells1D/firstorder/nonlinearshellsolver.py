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

def remove_fixed_nodes_vector(v, fixed_nodes_indicies, all_nodes_count, bc):
    indicies_to_exclude = i_exclude(fixed_nodes_indicies, all_nodes_count, bc)

    free_nodes1 = [i for i in range(v.shape[0]) if i not in indicies_to_exclude]
    return v[free_nodes1]


def extend_with_fixed_nodes(eig_vectors, fixed_nodes_indicies, all_nodes_count, bc):
    indicies_to_exclude = i_exclude(fixed_nodes_indicies, all_nodes_count, bc)
    res = eig_vectors
    for i in indicies_to_exclude:
        res = np.insert(res, i, 0, axis=0)

    return res


def i_exclude(fixed_nodes_indicies, nodes_count, bound_cond):
    fixed_indicies1 = [3 * x for x in fixed_nodes_indicies]
    fixed_indicies2 = []
    fixed_indicies3 = [3 * x + 2 for x in fixed_nodes_indicies]
    
    if (bound_cond == mod.FIXED_LEFT_RIGHT_EDGE):
        fixed_indicies1 = [3 * x for x in fixed_nodes_indicies]
        fixed_indicies2 = [3 * x + 1 for x in fixed_nodes_indicies]
        fixed_indicies3 = [3 * x + 2 for x in fixed_nodes_indicies]
    
    return sorted(fixed_indicies1+
                  fixed_indicies2+
                  fixed_indicies3)
    
    

def solve_nl(model, mesh, s_matrix, m_matrix, s_matrix_nl_1, s_matrix_nl_2, u_max, u_index = 0):

    s_ex = integrate_matrix(model, mesh, s_matrix)
    m_ex = integrate_matrix(model, mesh, m_matrix)
    

    fixed_nodes_indicies = mesh.get_fixed_nodes_indicies()

    s = remove_fixed_nodes(s_ex, fixed_nodes_indicies, mesh.nodes_count(), model.boundary_conditions)
    m = remove_fixed_nodes(m_ex, fixed_nodes_indicies, mesh.nodes_count(), model.boundary_conditions)


    lam, vec = la.eigh(s, m)

    vec_ex = extend_with_fixed_nodes(vec, fixed_nodes_indicies, mesh.nodes_count(), model.boundary_conditions)

    q = vec_ex[:,u_index]
    n = np.linalg.norm(q)
    
    w_max = get_max_w(q, mesh)
    
    res = normalize_w_only(q, u_max, w_max)
    
    s_nl_2_in = integrate_matrix_with_disp(model, mesh, s_matrix_nl_2, res)
    s_nl_2 = remove_fixed_nodes(s_nl_2_in, fixed_nodes_indicies, mesh.nodes_count(), model.boundary_conditions)
    
    s_nl_1_in = integrate_matrix_with_disp(model, mesh, s_matrix_nl_1, res)
    
    
    K = s + 0.75*s_nl_2
    
#    h = 0
#    for l in model.layers:
#        h = l.height()
#    
#    
#    A = u_max / h
#    print('=====koef 1D1O =====')
#    k1 = q.T.dot(s_ex).dot(q)
##    k2 = 0
##    if (A != 0):
#    k2 = q.T.dot(s_nl_2_in).dot(q)# / (A*A)
#    koef = (k1 + 0.75*k2)#/ lam[u_index])
#    print(vec[:,u_index].T.dot(K).dot(vec[:,u_index]))
#    print(koef)
#    print(vec[:,u_index].T.dot(m).dot(vec[:,u_index]))
#    
#    

    lam_nl = vec[:,u_index].T.dot(K).dot(vec[:,u_index])
    

#    
#    print('=====lam 1D1O =====')
#    
#    print(lam_nl)
    
#    lam2, vec2 = la.eigh(K, m)
#    
#    lam_nl = lam2[u_index]
#    
#    vec2 = extend_with_fixed_nodes(vec2, fixed_nodes_indicies, mesh.nodes_count(), model.boundary_conditions)
#
#    res2 = vec2[:,u_index]
#    
#    w_max2 = get_max_w(res2, mesh)
#    
#    res2 = normalize_w_only(res2, u_max, w_max2)
#    
#    print(la.norm(res-res2))
    
    
    b1 = -0.5*s_nl_1_in.dot(res)
    b2 = -0.25*s_nl_2_in.dot(res)
    
    b1 = remove_fixed_nodes_vector(b1, fixed_nodes_indicies, mesh.nodes_count(), model.boundary_conditions)
    b2 = remove_fixed_nodes_vector(b2, fixed_nodes_indicies, mesh.nodes_count(), model.boundary_conditions)
    
    U1 = la.solve(K, b1)
    U2 = la.solve(K - 4*lam_nl*m, b1)
    U3 = la.solve(K - 9*lam_nl*m, b2)
    
    U1 = extend_with_fixed_nodes(U1, fixed_nodes_indicies, mesh.nodes_count(), model.boundary_conditions)
    U2 = extend_with_fixed_nodes(U2, fixed_nodes_indicies, mesh.nodes_count(), model.boundary_conditions)
    U3 = extend_with_fixed_nodes(U3, fixed_nodes_indicies, mesh.nodes_count(), model.boundary_conditions)
    
    return lam_nl, res, U1, U2, U3, n
    

def integrate_matrix(model, mesh, matrix_func):
    N = 3 * (mesh.nodes_count())
    global_matrix = np.zeros((N, N))
    for element in mesh.elements:
        element_matrix = quadgch5nodes1dim(element_func, element, model.geometry, matrix_func)

        global_matrix += convertToGlobalMatrix(element_matrix, element, N)

    return global_matrix

def integrate_matrix_with_disp(model, mesh, matrix_func, disp):
    N = 3 * (mesh.nodes_count())
    global_matrix = np.zeros((N, N))
    for element in mesh.elements:
        u_element = matrices.get_u_element(element, disp, mesh.nodes_count())
        element_matrix = quadgch5nodes1dim(element_func_disp, element, model.geometry, matrix_func, u_element)

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


def element_func_disp(ksi, element, geometry, matrix_func, u_element):
    x1 = element.to_model_coordinates(ksi)
    x3 = 0
    x2 = 0
    
    u = matrices.get_u_deriv(element, u_element, x1, x2, x3)
    EM = matrix_func(element.material, geometry, x1, element.thickness, u)
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


def normalize(v, u_max):
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v*u_max / norm

def normalize_w_only(v, u_max, u3_max):
    if u3_max == 0:
        return v
    return v*u_max / u3_max

def get_max_w(v, mesh):
    w = v[range(2,3 * mesh.nodes_count(),3)]
    Wni = np.argmax(np.absolute(w))
    return w[Wni]
