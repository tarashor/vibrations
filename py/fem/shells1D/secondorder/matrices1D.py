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

def teta0(geometry, x1, x2, x3, h):
    A, K = geometry.get_A_and_K(x1, x2, x3)
    
    
    T = np.zeros((3, 3))
    
    Kh = K*h
    K2h2 = K*K*h*h
    
    T[0,0] = 16+6*Kh+K2h2
    T[0,1] = T[1,0] = 4*Kh+2*K2h2
    T[0,2] = T[2,0] = 16+16*Kh+4*K2h2
    T[1,1] = -2*Kh+K2h2
    T[0,2] = T[2,0] = -16*Kh+4*K2h2
    T[2,2] = -16+8*Kh+4*K2h2
    
    return T/32

def teta1(geometry, x1, x2, x3, h):
    A, K = geometry.get_A_and_K(x1, x2, x3)
    
    T = np.zeros((3, 3))
    
    Kh = K*h
    K2h2 = K*K*h*h
    
    T[0,0] = 2*Kh+K2h2
    T[0,1] = T[1,0] = -4*Kh+2*K2h2
    T[0,2] = T[2,0] = -16+4*K2h2
    T[1,1] = 16-6*Kh+K2h2
    T[0,2] = T[2,0] = 16-16*Kh+4*K2h2
    T[2,2] = -16-8*Kh+4*K2h2
    
    return T/32

def teta2(geometry, x1, x2, x3, h):
    A, K = geometry.get_A_and_K(x1, x2, x3)
    
    
    T = np.zeros((3, 3))
    
    Kh = K*h
    K2h2 = K*K*h*h
    
    T[0,0] = -4-4*Kh-K2h2
    T[0,1] = T[1,0] = 8-2*K2h2
    T[0,2] = T[2,0] = 16-8*Kh-K2h2
    T[1,1] = -4+4*Kh-K2h2
    T[0,2] = T[2,0] = 16+8*Kh-4*K2h2
    T[2,2] = 32-4*K2h2
    
    return T/32

def rotations_to_strain_nl(geometry, x1, x2, x3, h, u):
    T0 = teta0(geometry, x1, x2, x3, h)
    T1 = teta1(geometry, x1, x2, x3, h)
    T2 = teta2(geometry, x1, x2, x3, h)
    
    T11 = np.concatenate((T0, T1), axis=0)
    T11 = np.concatenate((T11, T2), axis=0)
    
    T33 = T11
    
    T13 = np.zeros((9, 3))
    
    T = np.concatenate((T11, T33), axis=0)
    T = np.concatenate((T, T13), axis=0)
    
    
    W = u_to_rotations(geometry,x1,0,0, h)
    
    rotations = W.dot(u)
    
    w = np.zeros((3, 9))
    
    w[0,0] = rotations[0]
    w[0,1] = rotations[1]
    w[0,2] = rotations[2]
    
    w[1,3] = rotations[0]
    w[1,4] = rotations[1]
    w[1,5] = rotations[2]
    
    w[2,6] = rotations[0]
    w[2,7] = rotations[1]
    w[2,8] = rotations[2]
    
    w_zero = np.zeros((3, 9))
    
    rot_l0 = np.concatenate((w, w_zero), axis=1)
    rot_l0 = np.concatenate((rot_l0, w_zero), axis=1)
    
    rot_l1 = np.concatenate((w_zero, w), axis=1)
    rot_l1 = np.concatenate((rot_l1, w_zero), axis=1)
    
    rot_l2 = np.zeros((3, 27))
    
    rot_l = np.concatenate((rot_l0, rot_l1), axis=0)
    rot_l = np.concatenate((rot_l, rot_l2), axis=0)
    

    return rot_l.dot(T).dot(W)



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
    u_nodes = np.zeros((12))
#    print(u[element.top_left_index])
    
    u10 = u[range(0,6 * nodes_count,6)]
    u11 = u[range(1,6 * nodes_count,6)]
    u12 = u[range(2,6 * nodes_count,6)]

    u30 = u[range(3,6 * nodes_count,6)]
    u31 = u[range(4,6 * nodes_count,6)]
    u32 = u[range(5,6 * nodes_count,6)]
    
    u_nodes[0] = u10[element.start_index]
    u_nodes[1] = u11[element.start_index]
    u_nodes[2] = u12[element.start_index]
    u_nodes[3] = u30[element.start_index]
    u_nodes[4] = u31[element.start_index]
    u_nodes[5] = u32[element.start_index]
        
    u_nodes[6] = u10[element.end_index]
    u_nodes[7] = u11[element.end_index]
    u_nodes[8] = u12[element.end_index]
    u_nodes[9] = u30[element.end_index]
    u_nodes[10] = u31[element.end_index]
    u_nodes[11] = u32[element.end_index]
    

    return u_nodes
    

def get_u_deriv(element, u_element, x1, x2, x3):
    h_e = element_aprox_functions(element, x1, x2, x3)

    return h_e.dot(u_element)


def get_C(material, geometry, x1, h):
    
    A, K = geometry.get_A_and_K(x1, 0, 0)
    
    C = material.matrix_C(geometry, x1, 0, 0)
    

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
    
    C11 = C[0,0]*h*C_
    C12 = C[0,2]*h*C_
    C21 = C[2,0]*h*C_
    C22 = C[2,2]*h*C_
    C33 = C[4,4]*h*C_
    
    C1 = C11
    C1 = np.concatenate((C1, C12), axis=1)
    C1 = np.concatenate((C1, np.zeros((3,3))), axis=1)
    
    
    C2 = C21
    C2 = np.concatenate((C2, C22), axis=1)
    C2 = np.concatenate((C2, np.zeros((3,3))), axis=1)
    
    C3 = np.zeros((3,6))
    C3 = np.concatenate((C3, C33), axis=1)
    
    C_M = np.concatenate((C1, C2), axis=0)
    C_M = np.concatenate((C_M, C3), axis=0)
    

    return C_M*A


def stiffness_matrix(material, geometry, x1, h):
    C = get_C(material, geometry, x1, h)
    
    E=u_to_strain(geometry,x1,0,0, h)
    
    
    return E.T.dot(C).dot(E)

def stiffness_matrix_nl_1(material, geometry, x1, h, u):
    C = get_C(material, geometry, x1, h)
    
    E=u_to_strain(geometry,x1,0,0, h)
    
    E_NL = rotations_to_strain_nl(geometry, x1, 0, 0, h, u)

    return 2*E_NL.T.dot(C).dot(E)+E.T.dot(C).dot(E_NL)

def stiffness_matrix_nl_2(material, geometry, x1, h, u):
    C = get_C(material, geometry, x1, h)
    
    E_NL = rotations_to_strain_nl(geometry, x1, 0, 0, h, u)

    return 2*E_NL.T.dot(C).dot(E_NL)
    


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




