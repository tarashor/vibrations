import numpy as np
from scipy import linalg as la
from scipy import integrate


def q2d(q, func, a1, b1, a2, b2, element, material, geometry, N):
    def ify(y):
        return q(func, a1, b1, args=(y, element, material, N))[0]
    return q(ify, a2, b2)[0]


def solve(model, mesh):
    s = stiffness_matrix(model, mesh)
    m = mass_matrix(model, mesh)

    # s = mesh.apply_boundary_conditions(s)
    # m = mesh.apply_boundary_conditions(m)
    lam, vec = la.eigh(s, m)
    return lam


def stiffness_matrix(model, mesh):
    N = 2 * (mesh.nodes_count())
    K = np.zeros((N, N))
    for element in mesh.elements:
        K_element = q2d(integrate.quad, k_element_func, -1, 1, -1, 1, element, mesh.material_for_element(element), model.geometry, N)
        K += K_element

    return K


def k_element_func(ksi, teta, element, material, geometry, N):
    k = np.zeros((N, N))
    C = const_matrix(material)
    E = grad_to_strain_linear_matrix()
    B = deriv_to_grad(alpha1, alpha2, geometry)
    I = ksiteta_to_alpha_matrix(element)
    H = lin_aprox_matrix(ksi, teta, element, N)
    J = jacobian(element)

    return k


def const_matrix(material):
    N = 6
    C = np.zeros((N, N))

    return C


def grad_to_strain_linear_matrix():
    E = np.zeros((6, 9))
    E[0, 0] = 1
    E[1, 4] = 1
    E[2, 8] = 1
    E[3, 1] = E[3, 3] = 1
    E[4, 2] = E[4, 6] = 1
    E[5, 5] = E[5, 7] = 1

    return E


def deriv_to_grad(alpha1, alpha2, geometry):
    B = np.zeros((9, 12))

    q = 1 + geometry.curvature * alpha2
    a = (np.pi + geometry.curvature * geometry.width) / 2 - geometry.curvature * alpha1
    w = q + geometry.corrugation_amplitude * geometry.curvature * np.cos(geometry.corrugation_frequency * a)
    z = geometry.corrugation_amplitude * geometry.corrugation_frequency * geometry.curvature * np.sin(geometry.corrugation_frequency * a)

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

    return B


def ksiteta_to_alpha_matrix(element):
    I = np.identity((6, 6))

    I[1, 1] = I[4, 4] = 2 / element.width()
    I[2, 2] = I[5, 5] = 2 / element.height()

    return I


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

    H = np.identity((6, N))
    d = N / 2
    H[0, element.top_left_index] = H[3, d + element.top_left_index] = f0

    return H


def mass_matrix(model, mesh):
    N = 2 * (mesh.nodes_count())
    M = np.zeros((N, N))
    return M

# b=np.array([[1,2], [3,4]])
# print(a.dot(b)==b)
