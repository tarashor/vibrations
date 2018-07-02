import numpy as np


def deriv_ksiteta_to_alpha(element):

    D = np.zeros((12, 12))

    for i in range(12):
        if (i % 2 == 0):
            D[i, i] = 1
        else:
            D[i, i] = 2 / element.width()

    return D


def lin_aprox_matrix(element, x1, x2, x3):

    ksi = element.to_element_coordinates(x1)

    f0 = 0.5 * (1 - ksi)
    f1 = 0.5 * (1 + ksi)

    Df0 = -0.5
    Df1 = 0.5

    H = np.zeros((12, 12))
    for i in range(6):
        H[2*i, i] = f0
        H[2*i, i+6] = f1
        H[2*i+1, i] = Df0
        H[2*i+1, i+6] = Df1
    
    return H


def element_aprox_functions(element, x1, x2, x3):
    H = deriv_ksiteta_to_alpha(element).dot(lin_aprox_matrix(element, x1, x2, x3))
    return H

def u_to_strain(geometry, x1, x2, x3, h):
    A, K = geometry.get_A_and_K(x1, x2, x3)
    
    hInv= 1/h
    
    E = np.zeros((9, 12))
    
    E[0, 1] = E[1, 3] = E[2, 5] = E[6, 7] = E[7, 9] = E[8, 11] = 1/A;
    E[0, 6] = E[1, 8] = E[2, 10] = E[8, 4] = K;
    E[5, 10] = 2 * K;

    E[3, 6] = -hInv + K / 2;
    E[3, 8] = E[6, 2] = E[7, 2] = hInv - K / 2;
    E[3, 10] = E[6, 4] = 4 * hInv - 2 * K;

    E[4, 6] = E[6, 0] = E[7, 0] = -hInv - K / 2;
    E[4, 8] = hInv + K / 2;
    E[4, 10] = E[7, 4] = -4 * hInv - 2 * K;

    return E

def u_to_rotations(geometry, x1, x2, x3, h):
    A, K = geometry.get_A_and_K(x1, x2, x3)
    
    hInv= 1/h
    
    W = np.zeros((3, 12))
    
    W[0, 0] = -hInv + 3*K / 2
    W[0, 2] = hInv - K / 2
    W[0, 4] = 4*hInv - 2*K
    
    W[0, 7] = -1/A
    
    W[1, 0] = -hInv - K / 2
    W[1, 2] = hInv + 3*K / 2
    W[1, 4] = -4*hInv - 2*K
    
    W[1, 9] = -1/A
    
    W[2, 4] = 3*K
    
    W[2, 11] = -1/A
    
    return 0.5*W

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


def get_u_element(element, u, nodes_count):
    u_nodes = np.zeros((8))
#    print(u[element.top_left_index])

    u_nodes[0] = u[element.top_left_index]
    u_nodes[1] = u[element.top_right_index]
    u_nodes[2] = u[element.bottom_right_index]
    u_nodes[3] = u[element.bottom_left_index]

    u_nodes[4] = u[element.top_left_index + nodes_count]
    u_nodes[5] = u[element.top_right_index + nodes_count]
    u_nodes[6] = u[element.bottom_right_index + nodes_count]
    u_nodes[7] = u[element.bottom_left_index + nodes_count]

    return u_nodes
    

def get_u_deriv(element,u_element, x1, x2, x3):
    h_e = element_aprox_functions(element, x1, x2, x3)

    return h_e.dot(u_element)

def get_grad_u(element,geometry,u_element, x1, x2, x3):
    B = deriv_to_grad(geometry, x1, x2, x3)
    h_e = element_aprox_functions(element, x1, x2, x3)

    return B.dot(h_e).dot(u_element)

def get_C(material, geometry, x1, h):
    
    C = material.matrix_C(geometry, x1, 0, 0)
    
    A = np.zeros((9,9))

    C_ = np.zeros((3,3))
    
    C_[0,0] = 1/3
    C_[0,1] = 1/6
    C_[0,2] = 1/3
    C_[1,0] = 1/6
    C_[1,1] = 1/3
    C_[1,2] = 1/3
    C_[2,0] = 1/3
    C_[2,1] = 1/3
    C_[2,2] = 8/15
    
    A11 = C[0,0]*h*C_
    A12 = C[0,2]*h*C_
    A21 = C[2,0]*h*C_
    A22 = C[2,2]*h*C_
    A33 = C[4,4]*h*C_
    
    A1 = A11
    A1 = np.concatenate((A1, A12), axis=1)
    A1 = np.concatenate((A1, np.zeros((3,3))), axis=1)
    
    
    A2 = A21
    A2 = np.concatenate((A2, A22), axis=1)
    A2 = np.concatenate((A2, np.zeros((3,3))), axis=1)
    
    A3 = np.zeros((3,6))
    A3 = np.concatenate((A3, A33), axis=1)
    
    A = np.concatenate((A1, A2), axis=0)
    A = np.concatenate((A, A3), axis=0)
    

    return A


def stiffness_matrix(material, geometry, x1, h):
    C = get_C(material, geometry, x1, h)
    
    
    E=u_to_strain(geometry,x1,0,0, h)
    
#    print(E.T)
    
#    print(C)
    
    return E.T.dot(C).dot(E)

def stiffness_matrix_nl_1(material, geometry, x1, h):
    E_NL_1 = deformations_nl_1(geometry, grad_u, x1, x2, x3)
    E_NL_2 = deformations_nl_2(geometry, grad_u, x1, x2, x3)
    C = material.matrix_C(geometry, x1, x2, x3)
    E = grad_to_strain()
    B = deriv_to_grad(geometry, x1, x2, x3)
    gj = geometry.getJacobian(x1, x2, x3)
    E_NL = E_NL_1+E_NL_2
#    print(grad_u)
    return (B.T.dot((E_NL).T).dot(C).dot(E).dot(B)+B.T.dot((E).T).dot(C).dot(E_NL_1).dot(B))* gj
#    return (B.T.dot((E_NL).T).dot(C).dot(E).dot(B))* gj
#    return (B.T.dot((E).T).dot(C).dot(E_NL_1).dot(B))* gj

def stiffness_matrix_nl_2(material, geometry, x1, x2, x3, grad_u):
    E_NL_1 = deformations_nl_1(geometry, grad_u, x1, x2, x3)
    E_NL_2 = deformations_nl_2(geometry, grad_u, x1, x2, x3)
    C = material.matrix_C(geometry, x1, x2, x3)
    B = deriv_to_grad(geometry, x1, x2, x3)
    gj = geometry.getJacobian(x1, x2, x3)
    E_NL = E_NL_1+E_NL_2
    
    return B.T.dot((E_NL).T).dot(C).dot(E_NL_1).dot(B)* gj

def stiffness_matrix_nl(material, geometry, x1, h):
    C = get_C(material, geometry, x1, h)
    
    
    E_NL=u_to_strain(geometry,x1,0,0, h)
    
#    print(E.T)
    
#    print(C)
    
    return E.T.dot(C).dot(E)

def mass_matrix2(material, geometry, x1, h):
    
    A, K = geometry.get_A_and_K(x1, 0, 0)
    rho = material.rho
    
    M=np.zeros((12,12))
    
    M[0,0] = 0.25*A*h*rho + 0.0833333333333333*h*(-1.0*A*K*h*rho + 1.0*A*rho)
    
    M[0,2] = 0.166666666666667*A*h*rho

    M[0,4] = 0.05*A*K*h**2*rho + 0.5*A*h*rho + 0.0833333333333333*h*(-1.0*A*K*h*rho - 2.0*A*rho)
    
    M[2,0] = 0.166666666666667*A*h*rho
    
    M[2,2] = 0.25*A*h*rho + 0.0833333333333333*h*(1.0*A*K*h*rho + 1.0*A*rho)
    
    M[2,4] = -0.05*A*K*h**2*rho + 0.5*A*h*rho + 0.0833333333333333*h*(1.0*A*K*h*rho - 2.0*A*rho)
    
    M[4,0] = 0.05*A*K*h**2*rho + 0.5*A*h*rho + 0.0833333333333333*h*(-1.0*A*K*h*rho - 2.0*A*rho)
    
    M[4,2] = -0.05*A*K*h**2*rho + 0.5*A*h*rho + 0.0833333333333333*h*(1.0*A*K*h*rho - 2.0*A*rho)
    
    M[4,4] = 8*A*h*rho/15
    
    M[6,6] = 0.25*A*h*rho + 0.0833333333333333*h*(-1.0*A*K*h*rho + 1.0*A*rho)
    
    M[6,8] = 0.166666666666667*A*h*rho
    
    M[6,10] = 0.05*A*K*h**2*rho + 0.5*A*h*rho + 0.0833333333333333*h*(-1.0*A*K*h*rho - 2.0*A*rho)
    
    M[8,6] = 0.166666666666667*A*h*rho
    
    M[8,8] = 0.25*A*h*rho + 0.0833333333333333*h*(1.0*A*K*h*rho + 1.0*A*rho)
    
    M[8,10] = -0.05*A*K*h**2*rho + 0.5*A*h*rho + 0.0833333333333333*h*(1.0*A*K*h*rho - 2.0*A*rho)
    
    M[10,6] = 0.05*A*K*h**2*rho + 0.5*A*h*rho + 0.0833333333333333*h*(-1.0*A*K*h*rho - 2.0*A*rho)
    
    M[10,8] = -0.05*A*K*h**2*rho + 0.5*A*h*rho + 0.0833333333333333*h*(1.0*A*K*h*rho - 2.0*A*rho)
    
    M[10,10] = 8*A*h*rho/15

    return M


def mass_matrix(material, geometry, x1, h):
    
    A, K = geometry.get_A_and_K(x1, 0, 0)
    rho = material.rho
    
    M=np.zeros((12,12))
    
    M[0,0] = 0.333333333333333*A*h*rho
    M[0,2] = 0.166666666666667*A*h*rho
    M[0,4] = 0.333333333333333*A*h*rho
    M[2,0] = 0.166666666666667*A*h*rho
    M[2,2] = 0.333333333333333*A*h*rho
    M[2,4] = 0.333333333333333*A*h*rho
    M[4,0] = 0.333333333333333*A*h*rho
    M[4,2] = 0.333333333333333*A*h*rho
    M[4,4] = 8*A*h*rho/15
    M[6,6] = 0.333333333333333*A*h*rho
    M[6,8] = 0.166666666666667*A*h*rho
    M[6,10] = 0.333333333333333*A*h*rho
    M[8,6] = 0.166666666666667*A*h*rho
    M[8,8] = 0.333333333333333*A*h*rho
    M[8,10] = 0.333333333333333*A*h*rho
    M[10,6] = 0.333333333333333*A*h*rho
    M[10,8] = 0.333333333333333*A*h*rho
    M[10,10] = 8*A*h*rho/15

    return M




