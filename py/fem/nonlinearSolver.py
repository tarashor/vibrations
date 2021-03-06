import numpy as np
# from . import mesh as m
from scipy import linalg as la


def solve(model, mesh):

    s = stiffness_matrix(model, mesh)
    m = mass_matrix(model, mesh)

    fixed_nodes_indicies = mesh.get_fixed_nodes_indicies()

    s = remove_fixed_nodes(s, fixed_nodes_indicies, mesh.nodes_count())
    m = remove_fixed_nodes(m, fixed_nodes_indicies, mesh.nodes_count())

    lam, vec = la.eigh(s, m)

    vec = extend_with_fixed_nodes(vec, fixed_nodes_indicies, mesh.nodes_count())

    return Result(lam, vec, mesh, model)


def solve_nonlinearity(model, mesh):

    s = stiffness_matrix(model, mesh)
    m = mass_matrix(model, mesh)

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


def stiffness_nl_matrix(res_prev, model, mesh):
    N = 2 * (mesh.nodes_count())
    K = np.zeros((N, N))
    for element in mesh.elements:
        material = mesh.material_for_element(element)
        K_element = quadgch5nodes2dim_prev_res(k_nl_element_func, element, material, model.geometry, res_prev)
        K += convertToGlobalMatrix(K_element, element, N)

    return K


def stiffness_matrix(model, mesh):
    N = 2 * (mesh.nodes_count())
    K = np.zeros((N, N))
    for element in mesh.elements:
        K_element = quadgch5nodes2dim(k_element_func, element, model.geometry)

        K += convertToGlobalMatrix(K_element, element, N)

    return K


def k_nl_element_func(ksi, teta, element, material, geometry, res):
    alpha1 = element.width() * ksi / 2 + (element.top_left.x + element.top_right.x) / 2
    alpha2 = element.height() * teta / 2 + (element.top_left.y + element.bottom_left.y) / 2
    C = const_matrix_isotropic(alpha1, alpha2, geometry, material)
    E = grad_to_strain_linear_matrix()
    B = deriv_to_grad(alpha1, alpha2, geometry)
    I_e = ksiteta_to_alpha_matrix(element)

    H = lin_aprox_matrix(ksi, teta)
    J = jacobian(element)

    u_e = np.zeros(8)

    for i in range(8):
        i_g = map_local_to_global_matrix_index(i, element, res.shape[0])
        u_e[i] = res[i_g]

    grad_u = B.dot(I_e).dot(H).dot(u_e)
    # print(grad_u)
    E_NL = grad_to_strain_nonlinear_matrix(alpha1, alpha2, geometry, grad_u)

    return H.T.dot(I_e.T).dot(B.T).dot(E_NL.T.dot(C).dot(E + E_NL) + E.T.dot(C).dot(E_NL)).dot(B).dot(I_e).dot(H) * J


def k_element_func(ksi, teta, element, material, geometry):
    alpha1 = element.width() * ksi / 2 + (element.top_left.x + element.top_right.x) / 2
    alpha2 = element.height() * teta / 2 + (element.top_left.y + element.bottom_left.y) / 2
    C = material.tensor_C(alpha1, alpha2, geometry, material)
    E = grad_to_strain_linear_matrix()
    B = deriv_to_grad(alpha1, alpha2, geometry)
    N = element_aprox_functions(element, alpha1, alpha2)
    J = jacobian(element)

    return N.T.dot(B.T).dot(E.T).dot(C).dot(E).dot(B).dot(N) * J


def grad_to_strain_linear_matrix():
    E = np.zeros((6, 9))
    E[0, 0] = 1
    E[1, 4] = 1
    E[2, 8] = 1
    E[3, 1] = E[3, 3] = 1
    E[4, 2] = E[4, 6] = 1
    E[5, 5] = E[5, 7] = 1

    return E


def grad_to_strain_nonlinear_matrix(alpha1, alpha2, geometry, grad_u):
    E = np.zeros((6, 9))
    g11 = geometry.get_g_11(alpha1, alpha2)

    d1u1 = grad_u[0]
    d2u1 = grad_u[1]
    d1u2 = grad_u[2]
    d2u2 = grad_u[3]

    E[0, 0] = g11 * d1u1
    E[0, 3] = d1u2
    E[1, 1] = g11 * d2u1
    E[1, 4] = d2u2

    E[3, 0] = g11 * d2u1
    E[3, 1] = g11 * d1u1
    E[3, 3] = d2u2
    E[3, 4] = d1u2

    E[4, 2] = g11 * d1u1
    E[4, 5] = d1u2

    E[5, 2] = g11 * d2u1
    E[5, 5] = d2u2

    return E


def deriv_to_grad(alpha1, alpha2, geometry):
    B = np.zeros((9, 12))

    G111, G211, G112, G121 = geometry.get_G(alpha1, alpha2)

    B[0, 0] = -G111
    B[0, 1] = 1
    B[0, 4] = -G211

    B[1, 0] = -G112
    B[1, 2] = 1

    B[2, 3] = 1

    B[3, 0] = -G121
    B[3, 5] = 1

    B[4, 6] = 1

    B[5, 7] = 1

    B[6, 9] = 1

    B[7, 10] = 1

    B[8, 11] = 1

    return B


def mass_matrix(model, mesh):
    N = 2 * (mesh.nodes_count())
    M = np.zeros((N, N))
    for element in mesh.elements:
        M_element = quadgch5nodes2dim(m_element_func, element, mesh.material_for_element(element), model.geometry)
        M += convertToGlobalMatrix(M_element, element, N)

    return M


def m_element_func(ksi, teta, element, material, geometry):
    alpha1 = element.width() * ksi / 2 + (element.top_left.x + element.top_right.x) / 2
    alpha2 = element.height() * teta / 2 + (element.top_left.y + element.bottom_left.y) / 2
    G = geometry.get_metric_matrix(alpha1, alpha2)
    B_s = deriv_to_vect(alpha1, alpha2, geometry)
    N = element_aprox_functions(element, alpha1, alpha2)
    J = jacobian(element)

    return material.rho * H.T.dot(I_e.T.dot(B_s.T.dot(G.dot(B_s.dot(I_e.dot(H)))))) * J


def deriv_to_vect(alpha1, alpha2, geometry):
    B = np.zeros((3, 12))

    B[0, 0] = B[1, 4] = B[2, 8] = 1

    return B


def quadgch5nodes2dim(f, element, material, geometry):
    order = 5
    w = [0.23692689, 0.47862867, 0.56888889, 0.47862867, 0.23692689]
    x = [-0.90617985, -0.53846931, 0, 0.53846931, 0.90617985]

    res = w[0] * w[0] * f(x[0], x[0], element, material, geometry)

    for i in range(order):
        for j in range(order):
            if (i != 0 or j != 0):
                res += w[i] * w[j] * f(x[i], x[j], element, material, geometry)

    return res


def quadgch5nodes2dim_prev_res(f, element, material, geometry, res_prev):
    order = 5
    w = [0.23692689, 0.47862867, 0.56888889, 0.47862867, 0.23692689]
    x = [-0.90617985, -0.53846931, 0, 0.53846931, 0.90617985]

    res = w[0] * w[0] * f(x[0], x[0], element, material, geometry, res_prev)

    for i in range(order):
        for j in range(order):
            if (i != 0 or j != 0):
                res += w[i] * w[j] * f(x[i], x[j], element, material, geometry, res_prev)

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
