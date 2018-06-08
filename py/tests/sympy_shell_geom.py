from sympy import *
from IPython.display import display
from sympy.vector import CoordSys3D
N = CoordSys3D('N')
x1, x2, x3 = symbols("x_1 x_2 x_3")
alpha1, alpha2, alpha3 = symbols("alpha_1 alpha_2 alpha_3")
R, L, ga, gv = symbols("R L g_a g_v")
init_printing()


#%%

a1 = pi / 2 + (L / 2 - alpha1)/R

x = R * cos(a1)
y = alpha2
z = R * sin(a1)

r = x*N.i + y*N.j + z*N.k

display(r)

#%%


r1 = r.diff(alpha1)
r2 = r.diff(alpha2)
k1 = trigsimp(r1.magnitude())
k2 = trigsimp(r2.magnitude())
r1 = r1/k1 
r2 = r2/k2
n = r1.cross(r2)
n = trigsimp(n.normalize())
n


# In[12]:

R_alpha=r+alpha3*n
R1 = R_alpha.diff(alpha1)
R1 = trigsimp(R1)
R2 = R_alpha.diff(alpha2)
R3 = R_alpha.diff(alpha3)


eps=trigsimp(R1.dot(R2.cross(R3)))
R_1=simplify(trigsimp(R2.cross(R3)/eps))
R_2=simplify(trigsimp(R3.cross(R1)/eps))
R_3=simplify(trigsimp(R1.cross(R2)/eps))

display(R1)
display(R2)
display(R3)


#%%

# Jacobi matrix:

dx1da1=R1.dot(N.i)
dx1da2=R2.dot(N.i)
dx1da3=R3.dot(N.i)

dx2da1=R1.dot(N.j)
dx2da2=R2.dot(N.j)
dx2da3=R3.dot(N.j)

dx3da1=R1.dot(N.k)
dx3da2=R2.dot(N.k)
dx3da3=R3.dot(N.k)

A=Matrix([[dx1da1, dx1da2, dx1da3], [dx2da1, dx2da2, dx2da3], [dx3da1, dx3da2, dx3da3]])
A = simplify(A)

A_inv = trigsimp(A**-1)
A_inv = simplify(trigsimp(A_inv))

#%%

# Metric tensor

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

display(G)
display(G_con)

#%%

DIM = 3

G_con_diff = MutableDenseNDimArray.zeros(DIM, DIM, DIM)
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            xdiff = alpha1
            if (k == 0): 
                xdiff = alpha1
            elif (k == 1):
                xdiff = alpha2
            elif (k == 2):
                xdiff = alpha3
            
            G_con_diff[i,j,k]=G_con[i,j].diff(xdiff)
            
display(G_con_diff)

#%%
DIM = 3



GK = MutableDenseNDimArray.zeros(DIM, DIM, DIM)
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            res = S(0)
            for m in range(DIM):
                res = res + G[m,k]*(G_con_diff[i,m,j]+G_con_diff[j,m,i]-G_con_diff[i,j,m])
            GK[i,j,k] = simplify(S(1)/S(2)*res)


display(GK)

