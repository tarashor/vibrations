import numpy as np

def deriv_ksiteta_to_alpha(element):

    D = np.zeros((6, 6))

    D[0, 0] = D[2, 2] = D[4, 4] = 1
    D[1, 1] = D[3, 3] = D[5, 5] = 2 / element.width()

    return D


def lin_aprox_matrix(element, x1, x2, x3):

    ksi = element.to_element_coordinates(x1)

    f0 = 0.5 * (1 - ksi)
    f1 = 0.5 * (1 + ksi)

    Df0 = -0.5
    Df1 = 0.5

    H = np.zeros((6, 6))
    H[0, 0] = H[2, 1] = H[4, 2] = f0
    H[1, 0] = H[3, 1] = H[5, 2] = Df0
    H[0, 3] = H[2, 4] = H[4, 5] = f1
    H[1, 3] = H[3, 4] = H[5, 5] = Df1
    
    return H


def element_aprox_functions(element, x1, x2, x3):
    H = deriv_ksiteta_to_alpha(element).dot(lin_aprox_matrix(element, x1, x2, x3))
    return H


def u_to_strain(geometry, x1, x2, x3):
    E = np.zeros((3,6))
    A, K = geometry.get_A_and_K(x1, x2, x3)
    
    E[0,1]=1/A
    
    E[0,4]=K
    
    E[1,3]=1/A

    E[2,0]=-K

    E[2,2]=1

    E[2,5]=1/A

    return E


def deriv_to_vect():
    B = np.zeros((3, 12))

    B[0, 0] = B[1, 4] = B[2, 8] = 1

    return B


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
    u1 = u[range(0,3 * nodes_count,3)]
    g = u[range(1,3 * nodes_count,3)]
    w = u[range(2,3 * nodes_count,3)]
    
    u_nodes = np.zeros((6))

    u_nodes[0] = u1[element.start_index]
    u_nodes[1] = g[element.start_index]
    u_nodes[2] = w[element.start_index]

    u_nodes[3] = u1[element.end_index]
    u_nodes[4] = g[element.end_index]
    u_nodes[5] = w[element.end_index]

    return u_nodes
    

def get_u_deriv(element,u_element, x1, x2, x3):
    h_e = element_aprox_functions(element, x1, x2, x3)

    return h_e.dot(u_element)

def get_grad_u(element,geometry,u_element, x1, x2, x3):
    B = deriv_to_grad(geometry, x1, x2, x3)
    h_e = element_aprox_functions(element, x1, x2, x3)

    return B.dot(h_e).dot(u_element)

def u_to_rotations(geometry, x1, x2, x3):
    A, K = geometry.get_A_and_K(x1, x2, x3)
    
    W = np.zeros((2, 6))
    
    W[0, 0] = K
    W[0, 2] = 1
    W[0, 5] = -1/A
    
    W[1, 2] = 2*K
    
    return 0.5*W

def teta0(geometry, x1, x2, x3):
    A, K = geometry.get_A_and_K(x1, x2, x3)
    
    T = np.zeros((2, 2))
    
    T[0,0] = 1
    
    return T/2

def teta1(geometry, x1, x2, x3):
    A, K = geometry.get_A_and_K(x1, x2, x3)
    
    T = np.zeros((2, 2))
    
    T[0,0] = -K
    T[0,1] = T[1,0] = 1
    
    return T/2


def rotations_to_strain_nl(geometry, x1, x2, x3, u):
    T11 = teta0(geometry, x1, x2, x3)
    
    Tk = teta1(geometry, x1, x2, x3)
    
    T13 = np.zeros((2, 2))
    
    T = np.concatenate((T11, Tk), axis=0)
    
    T = np.concatenate((T, T13), axis=0)
    
    W = u_to_rotations(geometry,x1,0,0)
    
    rotations = W.dot(u)
    
#    print('rotations============')
#    print(rotations)
    
    w = np.zeros((3, 6))
    
    w[0,0] = rotations[0]
    w[0,1] = rotations[1]
    
    w[1,2] = rotations[0]
    w[1,3] = rotations[1]
    
    return w.dot(T).dot(W)


def get_C(material, geometry, x1, h):
    
    A, K = geometry.get_A_and_K(x1, 0, 0)
    
    C = material.matrix_C(geometry, x1, 0, 0)
    
    C_ = np.zeros((3,3))
    
    C_[0,0] = h*C[0,0]
    C_[1,1] = (h**3)/12*C[0,0]  
    C_[2,2] = 5*h*C[4,4]/6
    

    return C_*A


def stiffness_matrix(material, geometry, x1, h):
    C = get_C(material, geometry, x1, h)
    
    
    E=u_to_strain(geometry,x1,0,0)
    
    return E.T.dot(C).dot(E)

def stiffness_matrix_nl_1(material, geometry, x1, h, u):
    C = get_C(material, geometry, x1, h)
    
    E=u_to_strain(geometry,x1,0,0)
    
    E_NL = rotations_to_strain_nl(geometry, x1, 0, 0, u)

    return 2*E_NL.T.dot(C).dot(E)+E.T.dot(C).dot(E_NL)

def stiffness_matrix_nl_2(material, geometry, x1, h, u):
    C = get_C(material, geometry, x1, h)
    
    E_NL = rotations_to_strain_nl(geometry, x1, 0, 0, u)

    return 2*E_NL.T.dot(C).dot(E_NL)
    


def mass_matrix(material, geometry, x1, h):
    
    A, K = geometry.get_A_and_K(x1, 0, 0)
    rho = material.rho
    h3=h**3
    
    M=np.zeros((6,6))
    M[0,0]=A*h*rho
    M[0,2]=A*K*h3*rho/12

    M[2,0]=A*K*h3*rho/12
    M[2,2]=A*h3*rho/12

    M[4,4]=A*h*rho
    
    return M

