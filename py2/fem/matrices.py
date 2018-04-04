import numpy as np


def deriv_ksiteta_to_alpha(element):

    D = np.zeros((12, 6))

    D[0, 0] = D[4, 3] = 1
    D[1, 1] = D[5, 4] = 2 / element.width()
    D[2, 2] = D[6, 5] = 2 / element.height()

    return D


def lin_aprox_matrix(element, x1, x2, x3):

    ksi, teta = element.to_element_coordinates(x1, x2)

    f0 = 0.25 * (1 - ksi) * (1 + teta)
    f1 = 0.25 * (1 + ksi) * (1 + teta)
    f2 = 0.25 * (1 + ksi) * (1 - teta)
    f3 = 0.25 * (1 - ksi) * (1 - teta)

    Df0_Dksi = -0.25 * (1 + teta)
    Df1_Dksi = 0.25 * (1 + teta)
    Df2_Dksi = 0.25 * (1 - teta)
    Df3_Dksi = -0.25 * (1 - teta)

    Df0_Dteta = 0.25 * (1 - ksi)
    Df1_Dteta = 0.25 * (1 + ksi)
    Df2_Dteta = -0.25 * (1 + ksi)
    Df3_Dteta = -0.25 * (1 - ksi)

    H = np.zeros((6, 8))
    H[0, 0] = H[3, 4] = f0
    H[0, 1] = H[3, 5] = f1
    H[0, 2] = H[3, 6] = f2
    H[0, 3] = H[3, 7] = f3

    H[1, 0] = H[4, 4] = Df0_Dksi
    H[1, 1] = H[4, 5] = Df1_Dksi
    H[1, 2] = H[4, 6] = Df2_Dksi
    H[1, 3] = H[4, 7] = Df3_Dksi

    H[2, 0] = H[5, 4] = Df0_Dteta
    H[2, 1] = H[5, 5] = Df1_Dteta
    H[2, 2] = H[5, 6] = Df2_Dteta
    H[2, 3] = H[5, 7] = Df3_Dteta

    return H


def element_aprox_functions(element, x1, x2, x3):
    H = deriv_ksiteta_to_alpha(element).dot(lin_aprox_matrix(element, x1, x2, x3))
    return H


def grad_to_strain():
    E = np.zeros((6, 9))
    E[0, 0] = 1
    E[1, 4] = 1
    E[2, 8] = 1
    E[3, 1] = E[3, 3] = 1
    E[4, 2] = E[4, 6] = 1
    E[5, 5] = E[5, 7] = 1

    return E


def deriv_to_grad(geometry, x1, x2, x3):
    B = np.zeros((9, 12))

    G = geometry.kristophel_symbols(x1, x2, x3)

    B[0, 0] = -G[0, 0, 0]
    B[0, 1] = 1
    B[0, 4] = -G[1, 0, 0]

    B[1, 0] = -G[0, 0, 1]
    B[1, 2] = 1
    B[1, 4] = -G[1, 0, 1]

    B[2, 3] = 1

    B[3, 0] = -G[0, 1, 0]
    B[3, 5] = 1
    B[3, 4] = -G[1, 1, 0]

    B[4, 6] = 1

    B[5, 7] = 1

    B[6, 9] = 1

    B[7, 10] = 1

    B[8, 11] = 1

#    K = geometry.curvature
#    q=1+K*x2
#
#    B[0, 0] = 0
#    B[0, 1] = 1/q
#    B[0, 4] = K/q
#
#    B[1, 2] = 1
#
#    B[2, 3] = 1
#
#    B[3, 0] = -K/q
#    B[3, 5] = 1/q
#
#    B[4, 6] = 1
#
#    B[5, 7] = 1
#
#    B[6, 9] = 1/q
#
#    B[7, 10] = 1
#
#    B[8, 11] = 1

    return B


def deriv_to_vect():
    B = np.zeros((3, 12))

    B[0, 0] = B[1, 4] = B[2, 8] = 1

    return B

def matrix_C(material, geometry, x1, x2, x3):
    N = 6

    C = np.zeros((N, N))

    lam = self.v * self.E / ((1 + self.v) * (1 - 2 * self.v))
    mu = self.E / ((1 + self.v) * 2)

    g = geometry.metric_tensor(x1, x2, x3)
    
    g = np.linalg.inv(g)

    for i in range(N):
        for j in range(N):
            n, m = self.__get_index_conv(i)
            k, l = self.__get_index_conv(j)
            C[i, j] = mu * (g[n, k] * g[m, l] + g[n, l] * g[m, k]) + lam * g[n, m] * g[k, l]

    return C

def get_index_conv(index):
    i = 0
    j = 0
    if (index == 0):
        i = 0
        j = 0
    elif (index == 1):
        i = 1
        j = 1
    elif (index == 2):
        i = 2
        j = 2
    elif (index == 3):
        i = 0
        j = 1
    elif (index == 4):
        i = 0
        j = 2
    elif (index == 5):
        i = 1
        j = 2

    return i, j

def model_stiffness_matrix(material, geometry, x1, x2, x3):
    K = zeros(12,12)
    
    return K

def model_mass_matrix(material, geometry, x1, x2, x3):
    M = zeros(12,12)
    
    return M
