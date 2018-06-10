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



def stiffness_matrix(material, geometry, x1, h):
    
    C = material.matrix_C(geometry, x1, 0, 0)
    
    C_ = np.zeros((3,3))
    
    C_[0,0] = h*C[0,0]
    C_[1,1] = (h**3)/12*C[0,0]
    C_[2,2] = h*C[4,4]
    
    E=u_to_strain(geometry,x1,0,0)
    
    return E.T.dot(C_).dot(E)


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

