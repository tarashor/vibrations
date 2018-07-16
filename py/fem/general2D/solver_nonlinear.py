from . import matrices2D as matrices
from . import result2D as r
import numpy as np
# from . import mesh as m
from scipy import linalg as la


def remove_fixed_nodes(matrix, fixed_nodes_indicies, all_nodes_count):
    indicies_to_exclude = i_exclude(fixed_nodes_indicies, all_nodes_count)

    free_nodes1 = [i for i in range(matrix.shape[0]) if i not in indicies_to_exclude]
    free_nodes2 = [i for i in range(matrix.shape[1]) if i not in indicies_to_exclude]
    return matrix[np.ix_(free_nodes1, free_nodes2)]

def remove_fixed_nodes_vector(v, fixed_nodes_indicies, all_nodes_count):
    indicies_to_exclude = i_exclude(fixed_nodes_indicies, all_nodes_count)

    free_nodes1 = [i for i in range(v.shape[0]) if i not in indicies_to_exclude]
    return v[free_nodes1]


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

def solve_nl(model, mesh, s_matrix, m_matrix, s_matrix_nl_1, s_matrix_nl_2, u_max, u_index = 0):

    s = integrate_matrix(model, mesh, s_matrix)
    m = integrate_matrix(model, mesh, m_matrix)
    

    fixed_nodes_indicies = mesh.get_fixed_nodes_indicies()

    s = remove_fixed_nodes(s, fixed_nodes_indicies, mesh.nodes_count())
    m = remove_fixed_nodes(m, fixed_nodes_indicies, mesh.nodes_count())


    lam, vec = la.eigh(s, m)

    vec_ex = extend_with_fixed_nodes(vec, fixed_nodes_indicies, mesh.nodes_count())

    q = vec_ex[:,u_index]
    u3_max = get_max_u3(q, mesh)
    
    res = normalize_u3_only(q, u_max, u3_max)
    
    s_nl_2_in = integrate_matrix_with_disp(model, mesh, s_matrix_nl_2, res)
    s_nl_2 = remove_fixed_nodes(s_nl_2_in, fixed_nodes_indicies, mesh.nodes_count())
    
    s_nl_1_in = integrate_matrix_with_disp(model, mesh, s_matrix_nl_1, res)
    
    K = s + 0.75*s_nl_2
    
#    h = 0
#    for l in model.layers:
#        h = l.height()
    
#    print('=====koef 2D =====')
#    koef = q.T.dot(s_nl_2_in).dot(q)
#    A = u_max / h
#    
#    print('A = {}'.format(A))
#    koef /= A*A
#    
#    print(koef / lam[u_index])
    
#    lam2, vec2 = la.eigh(K, m)
#    
#    
#    lam_nl = lam2[u_index]
#    
#    vec2 = extend_with_fixed_nodes(vec2, fixed_nodes_indicies, mesh.nodes_count())
#
#    res2 = vec2[:,u_index]
#    
#    u3_max2 = get_max_u3(res2, mesh)
#    
#    res2 = normalize_u3_only(res2, u_max, u3_max2)
    
#    print(la.norm(res-res2))
    
    lam_nl = vec[:,u_index].T.dot(K).dot(vec[:,u_index])
    
    k1 = vec_ex[:,u_index].T.dot(s_nl_1_in).dot(res)/lam_nl
    x1 = -0.5*k1
    x2 = k1/6
    k2=vec_ex[:,u_index].T.dot(s_nl_2_in).dot(res)/lam_nl
    x3 = k2/32
    print('x1={}'.format(x1*u3_max))
    print('x2={}'.format(x2*u3_max))
    print('x3={}'.format(x3*u3_max))
    
    
    b1 = -0.5*s_nl_1_in.dot(res)
    b2 = -0.25*s_nl_2_in.dot(res)
    
    b1 = remove_fixed_nodes_vector(b1, fixed_nodes_indicies, mesh.nodes_count())
    b2 = remove_fixed_nodes_vector(b2, fixed_nodes_indicies, mesh.nodes_count())
    
    U1 = la.solve(K, b1)
    U2 = la.solve(K - 4*lam_nl*m, b1)
    U3 = la.solve(K - 9*lam_nl*m, b2)
    
    U1 = extend_with_fixed_nodes(U1, fixed_nodes_indicies, mesh.nodes_count())
    U2 = extend_with_fixed_nodes(U2, fixed_nodes_indicies, mesh.nodes_count())
    U3 = extend_with_fixed_nodes(U3, fixed_nodes_indicies, mesh.nodes_count())
    
    return lam_nl, res, U1, U2, U3, u3_max


def convert_to_results(eigenvalues, eigenvectors, mesh, geometry):
    
    results = []
    for i in range(eigenvalues.size):
        freq = np.sqrt(eigenvalues[i])
        u1 = eigenvectors[:, i][0:mesh.nodes_count()]
        u3 = eigenvectors[:, i][mesh.nodes_count():2 * mesh.nodes_count()]
        u2 = np.zeros((mesh.nodes_count()))
        res = r.Result(freq, u1, u2, u3, mesh, geometry)
        results.append(res)

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
        u_element = matrices.get_u_element(element, disp, mesh.nodes_count())
        element_matrix = quadgch5nodes2dim(element_func_disp, element, model.geometry, matrix_func, u_element)

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


def normalize(v, u_max):
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v*u_max / norm

def normalize_u3_only(v, u_max, u3_max):
    if u3_max == 0:
        return v
    return v*u_max / u3_max

def get_max_u3(v, mesh):
    u3 = v[mesh.nodes_count():2 * mesh.nodes_count()]
    Wni = np.argmax(np.absolute(u3))
    print(Wni)
    return u3[Wni]

