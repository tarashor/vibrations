import numpy as np
from scipy import linalg as la


class Result:
    def __init__(self, lam, vec, mesh, model):
        self.lam = lam
        self.vec = vec
        self.mesh = mesh
        self.model = model

    def get_results_count(self):
        return self.lam

    def get_result(self, i):
        return self.lam[i], self.vec[:, i]


def solve(model, mesh):
    s = stiffness_matrix(model, mesh)
    m = mass_matrix(model, mesh)

    fixed_nodes_indicies = mesh.get_fixed_nodes_indicies(model.boundary_conditions)
    s = apply_boundary_conditions(s, fixed_nodes_indicies)
    m = apply_boundary_conditions(m, fixed_nodes_indicies)

    lam, vec = la.eigh(s, m)

    return Result(lam, vec, mesh)


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
        K_element = quadgch5nodes2dim(k_element_func, element, material, model.geometry, N)
        K += K_element

    return K


def k_element_func(ksi, teta, element, material, geometry, N):
    alpha1 = element.width() * ksi / 2 + (element.top_left.x + element.top_right.x) / 2
    alpha2 = element.height() * teta / 2 + (element.top_left.y + element.bottom_left.y) / 2
    C = const_matrix_isotropic(alpha1, alpha2, geometry, material)
    E = grad_to_strain_linear_matrix()
    B = deriv_to_grad(alpha1, alpha2, geometry)
    I_e = ksiteta_to_alpha_matrix(element)
    H = lin_aprox_matrix(ksi, teta, element, N)
    J = jacobian(element)

    return H.T.dot(I_e.T.dot(B.T.dot(E.T.dot(C.dot(E.dot(B.dot(I_e.dot(H)))))))) * J


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


def lin_aprox_matrix(ksi, teta, element, N):
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

    H = np.zeros((6, N))
    d = N // 2
    H[0, element.top_left_index] = f0
    H[3, d + element.top_left_index] = f0
    H[0, element.top_right_index] = H[3, d + element.top_right_index] = f1
    H[0, element.bottom_right_index] = H[3, d + element.bottom_right_index] = f2
    H[0, element.bottom_left_index] = H[3, d + element.bottom_left_index] = f3

    H[1, element.top_left_index] = H[4, d + element.top_left_index] = f0_ksi
    H[1, element.top_right_index] = H[4, d + element.top_right_index] = f1_ksi
    H[1, element.bottom_right_index] = H[4, d + element.bottom_right_index] = f2_ksi
    H[1, element.bottom_left_index] = H[4, d + element.bottom_left_index] = f3_ksi

    H[2, element.top_left_index] = H[5, d + element.top_left_index] = f0_teta
    H[2, element.top_right_index] = H[5, d + element.top_right_index] = f1_teta
    H[2, element.bottom_right_index] = H[5, d + element.bottom_right_index] = f2_teta
    H[2, element.bottom_left_index] = H[5, d + element.bottom_left_index] = f3_teta

    return H


def jacobian(element):
    return 0.25 * element.width() * element.height()


def mass_matrix(model, mesh):
    N = 2 * (mesh.nodes_count())
    M = np.zeros((N, N))
    for element in mesh.elements:
        M_element = quadgch5nodes2dim(m_element_func, element, mesh.material_for_element(element), model.geometry, N)
        M += M_element

    return M


def m_element_func(ksi, teta, element, material, geometry, N):
    alpha1 = element.width() * ksi / 2 + (element.top_left.x + element.top_right.x) / 2
    alpha2 = element.height() * teta / 2 + (element.top_left.y + element.bottom_left.y) / 2
    G = metric_matrix(alpha1, alpha2, geometry)
    B_s = deriv_to_vect(alpha1, alpha2, geometry)
    I_e = ksiteta_to_alpha_matrix(element)
    H = lin_aprox_matrix(ksi, teta, element, N)
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


def quadgch5nodes2dim(f, element, material, geometry, N):
    order = 5
    w = [0.23692689, 0.47862867, 0.56888889, 0.47862867, 0.23692689]
    x = [-0.90617985, -0.53846931, 0, 0.53846931, 0.90617985]

    res = np.zeros((N, N))

    for i in range(order):
        for j in range(order):
            res += w[i] * w[j] * f(x[i], x[j], element, material, geometry, N)

    return res

# b=np.array([[1,2], [3,4]])
# print(a.dot(b)==b)
