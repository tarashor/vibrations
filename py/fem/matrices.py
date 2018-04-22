import numpy as np


def deriv_ksiteta_to_alpha(element):

    D = np.zeros((12, 6))

    D[0, 0] = D[8, 3] = 1
    D[1, 1] = D[9, 4] = 2 / element.width()
    D[3, 2] = D[11, 5] = 2 / element.height()

    return D


def lin_aprox_matrix(element, x1, x2, x3):

    ksi, teta = element.to_element_coordinates(x1, x3)

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
    B[0, 8] = -G[1, 0, 0]

    B[1, 2] = 1
    
    B[2, 0] = -G[0, 0, 1]
    B[2, 3] = 1
    B[2, 8] = -G[1, 0, 1]

    B[3, 5] = 1

    B[4, 6] = 1

    B[5, 7] = 1

    B[6, 0] = -G[0, 1, 0]
    B[6, 9] = 1
    B[6, 8] = -G[1, 1, 0]

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

def tensor_C(material, geometry, x1, x2, x3):
    N = 6

    C = np.zeros((N, N))

    lam = material.lam()
    mu = material.mu()

    g = geometry.metric_tensor_inv(x1, x2, x3)

    for i in range(N):
        for j in range(N):
            n, m = get_index_conv(i)
            k, l = get_index_conv(j)
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

def deformations_nl(geometry, grad_u, x1, x2, x3):
    N = 3

    du = np.zeros((N, N))

    g = geometry.metric_tensor_inv(x1, x2, x3)

    for i in range(N):
        for j in range(N):
            index = i*N+j
            du[j,i] = grad_u[index]
    
    a_values = 0.5*du.dot(g)
    
    
    E_NL = np.zeros((6,9))
    E_NL[0,0] = a_values[0,0]
    E_NL[0,3] = a_values[0,1]
    E_NL[0,6] = a_values[0,2]
    
    E_NL[1,1] = a_values[1,0]
    E_NL[1,4] = a_values[1,1]
    E_NL[1,7] = a_values[1,2]
    
    E_NL[2,2] = a_values[2,0]
    E_NL[2,5] = a_values[2,1]
    E_NL[2,8] = a_values[2,2]
    
    E_NL[3,1] = 2*a_values[0,0]
    E_NL[3,4] = 2*a_values[0,1]
    E_NL[3,7] = 2*a_values[0,2]
    
    E_NL[4,2] = 2*a_values[0,0]
    E_NL[4,5] = 2*a_values[0,1]
    E_NL[4,8] = 2*a_values[0,2]
    
    E_NL[5,2] = 2*a_values[1,0]
    E_NL[5,5] = 2*a_values[1,1]
    E_NL[5,8] = 2*a_values[1,2]


    return E_NL



def stiffness_matrix(material, geometry, x1, x2, x3):
    C = tensor_C(material, geometry, x1, x2, x3)
    E = grad_to_strain()
    B = deriv_to_grad(geometry, x1, x2, x3)
    gj = geometry.getJacobian(x1, x2, x3)
    
    return B.T.dot(E.T).dot(C).dot(E).dot(B)* gj

def stiffness_matrix_nl(material, geometry, x1, x2, x3, grad_u):
    E_NL = deformations_nl(geometry, grad_u, x1, x2, x3)
    C = tensor_C(material, geometry, x1, x2, x3)
    E = grad_to_strain()
    B = deriv_to_grad(geometry, x1, x2, x3)
    gj = geometry.getJacobian(x1, x2, x3)
    
    return B.T.dot(E_NL.T).dot(C).dot(E).dot(B)* gj + B.T.dot(E.T).dot(C).dot(E_NL).dot(B)* gj

def mass_matrix(material, geometry, x1, x2, x3):
    g = geometry.metric_tensor(x1, x2, x3)
    g_inv = np.linalg.inv(g)
    gj = geometry.getJacobian(x1, x2, x3)
    
    B_s = deriv_to_vect()
    return material.rho * B_s.T.dot(g_inv.dot(B_s)) * gj
