from sympy import *
from sympy import MutableDenseNDimArray

DIM = 3

def getMetricTensorDownLame(H1, H2, H3):
    H12 = H1*H1
    H22 = H2*H2
    H32 = H3*H3
    
    G_con=Matrix([[H12, 0, 0],[0, H22, 0], [0, 0, H32]])
#    G_con=trigsimp(G_con)
    return G_con
    
    
def getMetricTensorUpLame(H1, H2, H3):
    
    H12 = H1*H1
    H22 = H2*H2
    H32 = H3*H3
    
    G=Matrix([[1/H12, 0, 0],[0, 1/H22, 0], [0, 0, 1/H32]])
#    G=trigsimp(G)
    return G

def getMetricTensorDown(R1, R2, R3):
    g_11=R1.dot(R1)
    g_12=R1.dot(R2)
    g_13=R1.dot(R3)
    
    g_21=R2.dot(R1)
    g_22=R2.dot(R2)
    g_23=R2.dot(R3)
    
    g_31=R3.dot(R1)
    g_32=R3.dot(R2)
    g_33=R3.dot(R3)
    
    G_con=Matrix([[g_11, g_12, g_13],[g_21, g_22, g_23], [g_31, g_32, g_33]])
    G_con=trigsimp(G_con)
    return G_con
    
    
def getMetricTensorUp(R_1, R_2, R_3):
    g11=R_1.dot(R_1)
    g12=R_1.dot(R_2)
    g13=R_1.dot(R_3)
    
    g21=R_2.dot(R_1)
    g22=R_2.dot(R_2)
    g23=R_2.dot(R_3)
    
    g31=R_3.dot(R_1)
    g32=R_3.dot(R_2)
    g33=R_3.dot(R_3)
    
    G=Matrix([[g11, g12, g13],[g21, g22, g23], [g31, g32, g33]])
    G=trigsimp(G)
    return G

def getChristoffelSymbols2(G_up, G_down_diff, axis):
    DIM = 3
    
    
    GK = MutableDenseNDimArray.zeros(DIM, DIM, DIM)
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                res = S(0)
                for m in range(DIM):
                    res = res + G_up[m,k]*(G_down_diff[i,m,j]+G_down_diff[j,m,i]-G_down_diff[i,j,m])
                GK[i,j,k] = S(1)/S(2)*res
    
    
    return GK

def getStiffnessTensor():
    
    C = MutableDenseNDimArray.zeros(DIM, DIM, DIM, DIM)
    
    for i in range(3):
        for j in range(3):        
            for k in range(3):
                for l in range(3):
                    elem_index = 'C^{{{}{}{}{}}}'.format(i+1, j+1, k+1, l+1)
                    el = Symbol(elem_index)
                    C[i,j,k,l] = el
                    
    return C

def getCIndecies(index):
    if (index == 0):
        return 0, 0
    elif (index == 1):
        return 1, 1
    elif (index == 2):
        return 2, 2
    elif (index == 3):
        return 0, 1
    elif (index == 4):
        return 0, 2
    elif (index == 5):
        return 1, 2

def getSymetricStiffnessTensor():
    
    C_x_symmetry = MutableDenseNDimArray.zeros(DIM, DIM, DIM, DIM)
    for s in range(6):
        for t in range(s, 6):
            i,j = getCIndecies(s)
            k,l = getCIndecies(t)
            elem_index = 'C^{{{}{}{}{}}}'.format(i+1, j+1, k+1, l+1)
            el = Symbol(elem_index)
            C_x_symmetry[i,j,k,l] = el
            C_x_symmetry[i,j,l,k] = el
            C_x_symmetry[j,i,k,l] = el
            C_x_symmetry[j,i,l,k] = el
            C_x_symmetry[k,l,i,j] = el
            C_x_symmetry[k,l,j,i] = el
            C_x_symmetry[l,k,i,j] = el
            C_x_symmetry[l,k,j,i] = el

                
    return C_x_symmetry

def convertStiffnessTensorToMatrix(C):

    C_matrix = zeros(6)
    
    for s in range(6):
        for t in range(6):
            i,j = getCIndecies(s)
            k,l = getCIndecies(t)
            C_matrix[s,t] = C[i,j,k,l]
            
    return C_matrix

def convertStiffnessMatrixToTensor(C):

    C_tensor = MutableDenseNDimArray.zeros(DIM, DIM, DIM, DIM)
            
    for s in range(6):
        for t in range(s, 6):
            i,j = getCIndecies(s)
            k,l = getCIndecies(t)
            el = C_tensor[s, t]
            C_tensor[i,j,k,l] = el
            C_tensor[i,j,l,k] = el
            C_tensor[j,i,k,l] = el
            C_tensor[j,i,l,k] = el
            C_tensor[k,l,i,j] = el
            C_tensor[k,l,j,i] = el
            C_tensor[l,k,i,j] = el
            C_tensor[l,k,j,i] = el
            
    return C_tensor

def getTensor4ElementInOtherCoordinates(C, A, q, p, s, t):
    res = S(0)
    for i in range(DIM):
        for j in range(DIM):        
            for k in range(DIM):
                for l in range(DIM):
                    res += C[i,j,k,l]*A[q,i]*A[p,j]*A[s,k]*A[t,l]
    return simplify(trigsimp(res))


def getTensor4InOtherCoordinates(C, A):
    C_alpha = MutableDenseNDimArray.zeros(DIM, DIM, DIM, DIM)

    for i in range(DIM):
        for j in range(DIM):        
            for k in range(DIM):
                for l in range(DIM):
                    c = getTensor4ElementInOtherCoordinates(C_isotropic, A, i, j, k, l)
                    C_alpha[i,j,k,l] = c
    
    return C_alpha

def getOrthotropicStiffnessTensor():

    C_orthotropic_matrix = zeros(6)
    
    for s in range(6):
        for t in range(s, 6):
            elem_index = 'C^{{{}{}}}'.format(s+1, t+1)
            el = Symbol(elem_index)
            if ((s < 3 and t < 3) or t == s):
                C_orthotropic_matrix[s,t] = el
                C_orthotropic_matrix[t,s] = el
            
    return convertStiffnessMatrixToTensor(C_orthotropic_matrix)

def getIsotropicStiffnessTensor():

    C_isotropic_matrix = zeros(6)

    mu = Symbol('mu')
    la = Symbol('lambda')
    
    for s in range(6):
        for t in range(s, 6):
            if (s < 3 and t < 3):
                if(t != s):
                    C_isotropic_matrix[s,t] = la
                    C_isotropic_matrix[t,s] = la
                else:
                    C_isotropic_matrix[s,t] = 2*mu+la
                    C_isotropic_matrix[t,s] = 2*mu+la
            elif (s == t):
                C_isotropic_matrix[s,t] = mu
                C_isotropic_matrix[t,s] = mu
            
    return convertStiffnessMatrixToTensor(C_isotropic_matrix)

def getUHat3D(alpha1,alpha2,alpha3):
    
    u1, u2, u3 = symbols("u_1, u_2, u_3")
    
    
    du = zeros(DIM,DIM)
    for i in range(DIM):
        for j in range(DIM):
            du[i,j]=Symbol('u_{{{},{}}}'.format(i+1,j+1))
    
    
    gu = zeros(12,1) 
    gu[0] = u1
    gu[1] = du[0,0]
    gu[2] = du[0,1]
    gu[3] = du[0,2]
    
    gu[4] = u2
    gu[5] = du[1,0]
    gu[6] = du[1,1]
    gu[7] = du[1,2]
    
    gu[8] = u3
    gu[9] = du[2,0]
    gu[10] = du[2,1]
    gu[11] = du[2,2]
    
    return gu

def getUHat3DPlane(alpha1,alpha2,alpha3):
    
    u1, u3 = symbols("u_1, u_3")
    u2=S(0)
    
    
    du = zeros(DIM,DIM)
    for i in range(DIM):
        for j in range(DIM):
            if (i != 1 and j != 1):
                du[i,j]=Symbol('u_{{{},{}}}'.format(i+1,j+1))
    
    
    gu = zeros(12,1) 
    gu[0] = u1
    gu[1] = du[0,0]
    gu[2] = du[0,1]
    gu[3] = du[0,2]
    
    gu[4] = u2
    gu[5] = du[1,0]
    gu[6] = du[1,1]
    gu[7] = du[1,2]
    
    gu[8] = u3
    gu[9] = du[2,0]
    gu[10] = du[2,1]
    gu[11] = du[2,2]
    
    return gu

def getUHatU3Main(alpha1,alpha2,alpha3):
    u1 = Function("u_1")
    u2 = Function("u_2")
    u3 = Function("u_3")
    
    
    gu = zeros(12,1) 
    gu[0] = u1(alpha1,alpha2,alpha3)
    
    gu[4] = u2(alpha1,alpha2,alpha3)
    
    gu[8] = u3(alpha1,alpha2,alpha3)
    gu[9] = u3(alpha1,alpha2,alpha3).diff(alpha1)
    gu[10] = u3(alpha1,alpha2,alpha3).diff(alpha2)
    
    return gu

#new
