import numpy as np

def ugw_to_u1u3(x1, x2, x3, h):
    p0=0.5-x3/h
    p1=0.5+x3/h
    p2=1-(2*x3/h)*(2*x3/h)
    
    dp0 = -1/h
    dp1 = 1/h
    dp2 = -8*x3/(h*h)
    
    L=np.zeros((12,12))
    
    L[0,0]=p0
    L[0,2]=p1
    L[0,4]=p2
    
    L[1,1]=p0
    L[1,3]=p1
    L[1,5]=p2
    
    L[3,0]=dp0
    L[3,2]=dp1
    L[3,4]=dp2
    
    L[8,6]=p0
    L[8,8]=p1
    L[8,10]=p2
    
    L[9,7]=p0
    L[9,9]=p1
    L[9,11]=p2
    
    L[11,6]=dp0
    L[11,8]=dp1
    L[11,10]=dp2

    return L


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

#    G[i, j, k]

    B[0, 0] = -G[0, 0, 0]
    B[0, 1] = 1
    B[0, 4] = -G[0, 0, 1]
    B[0, 8] = -G[0, 0, 2]
    
    B[1, 0] = -G[0, 1, 0]
    B[1, 2] = 1
    B[1, 4] = -G[0, 1, 1]
    B[1, 8] = -G[0, 1, 2]
    
    B[2, 0] = -G[0, 2, 0]
    B[2, 3] = 1
    B[2, 4] = -G[0, 2, 1]
    B[2, 8] = -G[0, 2, 2]
    
    B[3, 0] = -G[1, 0, 0]
    B[3, 5] = 1
    B[3, 4] = -G[1, 0, 1]
    B[3, 8] = -G[1, 0, 2]
    
    B[4, 0] = -G[1, 1, 0]
    B[4, 6] = 1
    B[4, 4] = -G[1, 1, 1]
    B[4, 8] = -G[1, 1, 2]
    
    B[5, 0] = -G[1, 2, 0]
    B[5, 7] = 1
    B[5, 4] = -G[1, 2, 1]
    B[5, 8] = -G[1, 2, 2]
    
    B[6, 0] = -G[2, 0, 0]
    B[6, 9] = 1
    B[6, 4] = -G[2, 0, 1]
    B[6, 8] = -G[2, 0, 2]
    
    B[7, 0] = -G[2, 1, 0]
    B[7, 10] = 1
    B[7, 4] = -G[2, 1, 1]
    B[7, 8] = -G[2, 1, 2]
    
    B[8, 0] = -G[2, 2, 0]
    B[8, 11] = 1
    B[8, 4] = -G[2, 2, 1]
    B[8, 8] = -G[2, 2, 2]

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


def deformations_nl_1(geometry, grad_u, x1, x2, x3):
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
    
    E_NL[4,0] = 2*a_values[2,0]
    E_NL[4,3] = 2*a_values[2,1]
    E_NL[4,6] = 2*a_values[2,2]
    
    E_NL[5,2] = 2*a_values[1,0]
    E_NL[5,5] = 2*a_values[1,1]
    E_NL[5,8] = 2*a_values[1,2]


    return E_NL

def deformations_nl_2(geometry, grad_u, x1, x2, x3):
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
    
    E_NL[3,0] = 2*a_values[1,0]
    E_NL[3,3] = 2*a_values[1,1]
    E_NL[3,6] = 2*a_values[1,2]
    
    E_NL[4,2] = 2*a_values[0,0]
    E_NL[4,5] = 2*a_values[0,1]
    E_NL[4,8] = 2*a_values[0,2]
    
    E_NL[5,1] = 2*a_values[2,0]
    E_NL[5,4] = 2*a_values[2,1]
    E_NL[5,7] = 2*a_values[2,2]


    return E_NL


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



def stiffness_matrix(material, geometry, x1, x2, x3):
    C = tensor_C(material, geometry, x1, x2, x3)
    E = grad_to_strain()
    B = deriv_to_grad(geometry, x1, x2, x3)
    gj = geometry.getJacobian(x1, x2, x3)
    
    return B.T.dot(E.T).dot(C).dot(E).dot(B)* gj

def stiffness_matrix_nl(material, geometry, x1, x2, x3, grad_u):
    E_NL_1 = deformations_nl_1(geometry, grad_u, x1, x2, x3)
    E_NL_2 = deformations_nl_2(geometry, grad_u, x1, x2, x3)
    C = tensor_C(material, geometry, x1, x2, x3)
    E = grad_to_strain()
    B = deriv_to_grad(geometry, x1, x2, x3)
    gj = geometry.getJacobian(x1, x2, x3)
    E_NL = E_NL_1+E_NL_2
    
    return B.T.dot((E_NL).T).dot(C).dot(E_NL_1).dot(B)* gj

def mass_matrix(material, geometry, x1, x2, x3):
    g = geometry.metric_tensor(x1, x2, x3)
    g_inv = np.linalg.inv(g)
    gj = geometry.getJacobian(x1, x2, x3)
    
    B_s = deriv_to_vect()
    return material.rho * B_s.T.dot(g_inv.dot(B_s)) * gj

