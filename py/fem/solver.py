import numpy as np
# from . import mesh as m
from scipy import linalg as la


class Result:
    def __init__(self, lam, vec, mesh, model):
        self.lam = lam
        self.vec = vec
        self.mesh = mesh
        self.model = model

    def get_results_count(self):
        return np.sqrt(self.lam)

    def get_result(self, i):
        return np.sqrt(self.lam[i]), self.vec[:, i]


def solve(model, mesh):
    print("==================Solver==================")
    print("STARTED:")
    
    fixed_nodes_indicies = mesh.get_fixed_nodes_indicies()
    print("Fixed nodes: {}".format(fixed_nodes_indicies))
    
    print("===Stiffness matrix: STARTED===")
    s = stiffness_matrix(model, mesh)
    print("===Stiffness matrix: FINISHED===")
    print("===MASS matrix: STARTED===")
    m = mass_matrix(model, mesh)
    print("===MASS matrix: STARTED===")

    s = apply_boundary_conditions(s, fixed_nodes_indicies)
    m = apply_boundary_conditions(m, fixed_nodes_indicies)

    lam, vec = la.eigh(s, m)
    print("FINISHED")
    print("==================Solver==================")
    return Result(lam, vec, mesh, model)


def apply_boundary_conditions(matrix, fixed_nodes_indicies):
    if (matrix.shape[0] == matrix.shape[1]):
        all_nodes_count = matrix.shape[0]
        free_nodes = [i for i in range(all_nodes_count) if i not in fixed_nodes_indicies]
        return matrix[np.ix_(free_nodes, free_nodes)]


def stiffness_matrix(model, mesh):
    N = 2 * (mesh.nodes_count())
    K = np.zeros((N, N))
    for element in mesh.elements:
        material = mesh.material_for_element(element)
        K_element = quadgch5nodes2dim(k_element_func, element, material, model.geometry)
        
        K += convertToGlobalMatrix(K_element, element, N)

    return K


def k_element_func(ksi, teta, element, material, geometry):
    alpha1 = element.width() * ksi / 2 + (element.top_left.x + element.top_right.x) / 2
    alpha2 = element.height() * teta / 2 + (element.top_left.y + element.bottom_left.y) / 2
    C = const_matrix_isotropic(alpha1, alpha2, geometry, material)
    E = grad_to_strain_linear_matrix()
    B = deriv_to_grad(alpha1, alpha2, geometry)
    I_e = ksiteta_to_alpha_matrix(element)
    
    SMALL_I = I_e.T.dot(B.T).dot(E.T).dot(C).dot(E).dot(B).dot(I_e)
        
    H = lin_aprox_matrix(ksi, teta)
    J = jacobian(element)

    return H.T.dot(SMALL_I).dot(H) * J


def const_matrix_isotropic(alpha1, alpha2, geometry, material):
    C = np.zeros((6, 6))
    g11 = get_g_11(alpha1, alpha2, geometry)
    v = material.v

    C[0, 0] = g11 * g11 * (1 - v)
    C[1, 1] = 1 - v
    C[2, 2] = 1 - v
    C[0, 1] = C[1, 0] = g11 * v
    C[0, 2] = C[2, 0] = g11 * v
    C[1, 2] = v
    C[2, 1] = v

    C[3, 3] = g11 * (1 - 2 * v) * 0.5
    C[4, 4] = g11 * (1 - 2 * v) * 0.5
    C[5, 5] = (1 - 2 * v) * 0.5

    koef = material.E / ((1 + v) * (1 - 2 * v))

    return koef * C[np.ix_([0, 1, 3], [0, 1, 3])]


def grad_to_strain_linear_matrix():
    E = np.zeros((6, 9))
    E[0, 0] = 1
    E[1, 4] = 1
    E[2, 8] = 1
    E[3, 1] = E[3, 3] = 1
    E[4, 2] = E[4, 6] = 1
    E[5, 5] = E[5, 7] = 1

    return E[np.ix_([0, 1, 3], [0, 1, 3, 4])]


def deriv_to_grad(alpha1, alpha2, geometry):
    B = np.zeros((9, 12))

    q, a, w, z = get_metric_tensor_components(alpha1, alpha2, geometry)

    G111 = 2 * z * geometry.curvature / w
    G211 = -geometry.corrugation_amplitude * geometry.corrugation_frequency * geometry.corrugation_frequency * geometry.curvature * geometry.curvature * np.cos(geometry.corrugation_frequency * a) - w * geometry.curvature - 2 * z * z * geometry.curvature / w
    G112 = geometry.curvature / q
    G121 = G112

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

    return B[np.ix_([0, 1, 3, 4], [0, 1, 2, 4, 5, 6])]


def ksiteta_to_alpha_matrix(element):
    I_e = np.identity(6)

    I_e[1, 1] = I_e[4, 4] = 2 / element.width()
    I_e[2, 2] = I_e[5, 5] = 2 / element.height()

    return I_e


def lin_aprox_matrix(ksi, teta):
    f0 = 0.25 * (1 - ksi) * (1 + teta)
    f1 = 0.25 * (1 + ksi) * (1 + teta)
    f2 = 0.25 * (1 + ksi) * (1 - teta)
    f3 = 0.25 * (1 - ksi) * (1 - teta)

    f0_ksi = -0.25 * (1 + teta)
    f1_ksi = 0.25 * (1 + teta)
    f2_ksi = 0.25 * (1 - teta)
    f3_ksi = -0.25 * (1 - teta)

    f0_teta = 0.25 * (1 - ksi)
    f1_teta = 0.25 * (1 + ksi)
    f2_teta = -0.25 * (1 + ksi)
    f3_teta = -0.25 * (1 - ksi)

    H = np.zeros((6, 8))
    H[0, 0] = H[3, 4] = f0
    H[0, 1] = H[3, 5] = f1
    H[0, 2] = H[3, 6] = f2
    H[0, 3] = H[3, 7] = f3

    H[1, 0] = H[4, 4] = f0_ksi
    H[1, 1] = H[4, 5] = f1_ksi
    H[1, 2] = H[4, 6] = f2_ksi
    H[1, 3] = H[4, 7] = f3_ksi

    H[2, 0] = H[5, 4] = f0_teta
    H[2, 1] = H[5, 5] = f1_teta
    H[2, 2] = H[5, 6] = f2_teta
    H[2, 3] = H[5, 7] = f3_teta

    return H


def jacobian(element):
    return 0.25 * element.width() * element.height()

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
    G = metric_matrix(alpha1, alpha2, geometry)
    B_s = deriv_to_vect(alpha1, alpha2, geometry)
    I_e = ksiteta_to_alpha_matrix(element)
    H = lin_aprox_matrix(ksi, teta)
    J = jacobian(element)

    return material.rho * H.T.dot(I_e.T.dot(B_s.T.dot(G.dot(B_s.dot(I_e.dot(H)))))) * J


def metric_matrix(alpha1, alpha2, geometry):

    G = np.zeros((3, 3))
    G[0, 0] = get_g_11(alpha1, alpha2, geometry)
    G[1, 1] = 1
    G[2, 2] = 1

    return G[np.ix_([0, 1], [0, 1])]


def deriv_to_vect(alpha1, alpha2, geometry):
    B = np.zeros((3, 12))

    B[0, 0] = B[1, 4] = B[2, 8] = 1

    return B[np.ix_([0, 1], [0, 1, 2, 4, 5, 6])]


def get_metric_tensor_components(alpha1, alpha2, geometry):
    q = 1 + geometry.curvature * alpha2
    a = (np.pi + geometry.curvature * geometry.width) / 2 - geometry.curvature * alpha1
    w = q + geometry.corrugation_amplitude * geometry.curvature * np.cos(geometry.corrugation_frequency * a)
    z = geometry.corrugation_amplitude * geometry.corrugation_frequency * geometry.curvature * np.sin(geometry.corrugation_frequency * a)
    return q, a, w, z


def get_g_11(alpha1, alpha2, geometry):
    q, a, w, z = get_metric_tensor_components(alpha1, alpha2, geometry)
    return 1 / (w * w + z * z)


def quadgch5nodes2dim(f, element, material, geometry):
    order = 5
    w = [0.23692689, 0.47862867, 0.56888889, 0.47862867, 0.23692689]
    x = [-0.90617985, -0.53846931, 0, 0.53846931, 0.90617985]

    res = w[0] * w[0] * f(x[0], x[0], element, material, geometry)

    for i in range(order):
        for j in range(order):
            if (i != 0 and j != 0):
                res += w[i] * w[j] * f(x[i], x[j], element, material, geometry)

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
    global_matrix = np.zeros((N,N))
    rows, columns = local_matrix.shape
    for i in range(rows):
        for j in range(columns):
            i_global = map_local_to_global_matrix_index(i, element, N)
            j_global = map_local_to_global_matrix_index(j, element, N)
            global_matrix[i_global, j_global] = local_matrix[i,j]
            
    return global_matrix

# b=np.array([[1,2], [3,4]])
# print(a.dot(b)==b)
