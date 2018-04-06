
# coding: utf-8

# # Shells

# ## Init symbols for *sympy*

# In[22]:


from sympy import *
from sympy.vector import CoordSys3D
import matplotlib.pyplot as plt
from matplotlib import cm
from IPython.display import display

N = CoordSys3D('N')
x1, x2, x3 = symbols("x_1 x_2 x_3")
alpha1, alpha2, alpha3 = symbols("alpha_1 alpha_2 alpha_3")
R, L, ga, gv = symbols("R L g_a g_v")
init_printing()


# ## Cylindrical coordinates

# In[23]:


a1 = pi / 2 + (L / 2 - alpha1)/R

x = R * cos(a1)
y = alpha2
z = R * sin(a1)

r = x*N.i + y*N.j + z*N.k

r1 = r.diff(alpha1)
r2 = r.diff(alpha2)
k1 = trigsimp(r1.magnitude())
k2 = trigsimp(r2.magnitude())
r1 = r1/k1 
r2 = r2/k2


n = r1.cross(r2)
n = trigsimp(n.normalize())


# ### Base Vectors $\vec{R}_1, \vec{R}_2, \vec{R}_3$

R_alpha=r+alpha3*n

R1=R_alpha.diff(alpha1)
R2=R_alpha.diff(alpha2)
R3=R_alpha.diff(alpha3)

display(trigsimp(R1))
display(R2)
display(R3)


# ### Draw

alpha_x = lambdify([R, L, alpha1, alpha3], R_alpha.dot(N.i), "numpy")
alpha_z = lambdify([R, L, alpha1, alpha3], R_alpha.dot(N.k), "numpy")

R_num = 1/0.8
L_num = 2

x1_start = 0
x1_end = L_num
x3_start = -0.05
x3_end = 0.05
plot_x1_elements = 100

dx1 = (x1_end - x1_start) / plot_x1_elements
X_init = []
Y_init = []
x2 = 0
x3 = 0

for i in range(plot_x1_elements + 1):
    x1 = x1_start + i * dx1

    x=alpha_x(R_num, L_num, x1, x3)
    z=alpha_z(R_num, L_num, x1, x3)

    X_init.append(x)
    Y_init.append(z)



plt.plot(X_init, Y_init, "r", label="початкова конфігурація")

geometry_title = "K={}".format(1/R_num)
plot_title = r"Форма панелі $L={}, h={}$".format(x1_end - x1_start, x3_end - x3_start)
if (len(geometry_title) > 0):
    plot_title = r"Форма панелі $L={}, h={}, {}$".format(x1_end - x1_start, x3_end - x3_start, geometry_title)

plt.title(plot_title)
plt.axes().set_aspect('equal', 'datalim')
plt.legend(loc='best')
plt.xlabel(r"$x_1$, м", fontsize=12)
plt.ylabel(r"$x_3$, м", fontsize=12)
plt.grid()
plt.show()



eps=trigsimp(R1.dot(R2.cross(R3)))
R_1=simplify(trigsimp(R2.cross(R3)/eps))
R_2=simplify(trigsimp(R3.cross(R1)/eps))
R_3=simplify(trigsimp(R1.cross(R2)/eps))


display(R_1)
display(R_2)
display(R_3)


#%%
# #### Jacobi matrix:

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
A=simplify(A)


# In[43]:


A_inv = trigsimp(A**-1)
A_inv = simplify(A_inv)


# In[44]:


jacobian = trigsimp(A.det())


# ### Metric tensor

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
G



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




# ### Derivatives of vectors



B = zeros(9, 12)
B[0,1] = (1+alpha3/R)**2
B[0,8] = (1+alpha3/R)/R

B[1,2] = (1+alpha3/R)**2

B[2,0] = (1+alpha3/R)/R
B[2,3] = (1+alpha3/R)**2

B[3,5] = S(1)
B[4,6] = S(1)
B[5,7] = S(1)

B[6,9] = S(1)
B[6,0] = -(1+alpha3/R)/R
B[7,10] = S(1)
B[8,11] = S(1)

#display(B)



B_con = zeros(9, 12)
B_con[0,1] = 1
B_con[0,8] = (1+alpha3/R)/R

B_con[1,2] = 1

B_con[2,0] = -1/(R*(1+alpha3/R))
B_con[2,3] = 1

B_con[3,5] = S(1)
B_con[4,6] = S(1)
B_con[5,7] = S(1)

B_con[6,0] = -1/(R*(1+alpha3/R))
B_con[6,9] = S(1)
B_con[7,10] = S(1)
B_con[8,11] = S(1)

koef=g_11.diff(alpha3)


# In[62]:


u_down_to_up=eye(12)
u_down_to_up[0,0]=g_11
u_down_to_up[1,1]=g_11
u_down_to_up[2,2]=g_11
u_down_to_up[3,3]=g_11
u_down_to_up[3,0]=koef


# In[63]:


#display(simplify(B_con*u_down_to_up))


# ### Deformations tensor

# In[64]:


E=zeros(6,9)
E[0,0]=1
E[1,4]=1
E[2,8]=1
E[3,1]=1
E[3,3]=1
E[4,2]=1
E[4,6]=1
E[5,5]=1
E[5,7]=1




# ### Elasticity tensor(stiffness tensor)
# 
# #### General form

# In[68]:


from sympy import MutableDenseNDimArray
C_x = MutableDenseNDimArray.zeros(3, 3, 3, 3)

for i in range(3):
    for j in range(3):        
        for k in range(3):
            for l in range(3):
                elem_index = 'C^{{{}{}{}{}}}'.format(i+1, j+1, k+1, l+1)
                el = Symbol(elem_index)
                C_x[i,j,k,l] = el
                
C_x


# #### Include symmetry

# In[69]:


C_x_symmetry = MutableDenseNDimArray.zeros(3, 3, 3, 3)

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

                
C_x_symmetry


# #### Isotropic material

# In[70]:


C_isotropic = MutableDenseNDimArray.zeros(3, 3, 3, 3)

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
            
for s in range(6):
    for t in range(s, 6):
        i,j = getCIndecies(s)
        k,l = getCIndecies(t)
        el = C_isotropic_matrix[s, t]
        C_isotropic[i,j,k,l] = el
        C_isotropic[i,j,l,k] = el
        C_isotropic[j,i,k,l] = el
        C_isotropic[j,i,l,k] = el
        C_isotropic[k,l,i,j] = el
        C_isotropic[k,l,j,i] = el
        C_isotropic[l,k,i,j] = el
        C_isotropic[l,k,j,i] = el

                
C_isotropic


# In[71]:


def getCalpha(C, A, q, p, s, t):
    res = S(0)
    for i in range(3):
        for j in range(3):        
            for k in range(3):
                for l in range(3):
                    res += C[i,j,k,l]*A[q,i]*A[p,j]*A[s,k]*A[t,l]
    return simplify(trigsimp(res))
                    


C_isotropic_alpha = MutableDenseNDimArray.zeros(3, 3, 3, 3)

for i in range(3):
    for j in range(3):        
        for k in range(3):
            for l in range(3):
                c = getCalpha(C_isotropic, A_inv, i, j, k, l)
                C_isotropic_alpha[i,j,k,l] = c

C_isotropic_alpha[0,0,0,0]


# In[72]:


C_isotropic_matrix_alpha = zeros(6)

for s in range(6):
    for t in range(6):
        i,j = getCIndecies(s)
        k,l = getCIndecies(t)
        C_isotropic_matrix_alpha[s,t] = C_isotropic_alpha[i,j,k,l]
        
C_isotropic_matrix_alpha


# #### Orthotropic material

# In[73]:


C_orthotropic = MutableDenseNDimArray.zeros(3, 3, 3, 3)

C_orthotropic_matrix = zeros(6)

for s in range(6):
    for t in range(s, 6):
        elem_index = 'C^{{{}{}}}'.format(s+1, t+1)
        el = Symbol(elem_index)
        if ((s < 3 and t < 3) or t == s):
            C_orthotropic_matrix[s,t] = el
            C_orthotropic_matrix[t,s] = el
            
for s in range(6):
    for t in range(s, 6):
        i,j = getCIndecies(s)
        k,l = getCIndecies(t)
        el = C_orthotropic_matrix[s, t]
        C_orthotropic[i,j,k,l] = el
        C_orthotropic[i,j,l,k] = el
        C_orthotropic[j,i,k,l] = el
        C_orthotropic[j,i,l,k] = el
        C_orthotropic[k,l,i,j] = el
        C_orthotropic[k,l,j,i] = el
        C_orthotropic[l,k,i,j] = el
        C_orthotropic[l,k,j,i] = el

                
C_orthotropic


# #### Orthotropic material in shell coordinates

                    


C_orthotropic_alpha = MutableDenseNDimArray.zeros(3, 3, 3, 3)

for i in range(3):
    for j in range(3):        
        for k in range(3):
            for l in range(3):
                c = getCalpha(C_orthotropic, A_inv, i, j, k, l)
                C_orthotropic_alpha[i,j,k,l] = c

C_orthotropic_alpha[0,0,0,0]


# In[75]:


C_orthotropic_matrix_alpha = zeros(6)

for s in range(6):
    for t in range(6):
        i,j = getCIndecies(s)
        k,l = getCIndecies(t)
        C_orthotropic_matrix_alpha[s,t] = C_orthotropic_alpha[i,j,k,l]
        
C_orthotropic_matrix_alpha


# ### Physical coordinates

# $u^1=\frac{u_{[1]}}{1+\frac{\alpha_3}{R}}$

# $\frac{\partial u^1} {\partial \alpha_3}=\frac{1}{1+\frac{\alpha_3}{R}} \frac{\partial u_{[1]}} {\partial \alpha_3} + u_{[1]} \frac{\partial} {\partial \alpha_3} \left( \frac{1}{1+\frac{\alpha_3}{R}} \right) = =\frac{1}{1+\frac{\alpha_3}{R}} \frac{\partial u_{[1]}} {\partial \alpha_3} - u_{[1]} \frac{1}{R \left( 1+\frac{\alpha_3}{R} \right)^2} $

# In[76]:


P=eye(12,12)
P[0,0]=1/(1+alpha3/R)
P[1,1]=1/(1+alpha3/R)
P[2,2]=1/(1+alpha3/R)
P[3,0]=-1/(R*(1+alpha3/R)**2)
P[3,3]=1/(1+alpha3/R)
P


# In[77]:


Def=simplify(E*B*P)
Def


# In[78]:


rows, cols = Def.shape
D_p=zeros(rows, cols)
q = 1+alpha3/R
for i in range(rows):
    ratio = 1
    if (i==0):
        ratio = q*q
    elif (i==3 or i == 4):
        ratio = q
    
    for j in range(cols):
        D_p[i,j] = Def[i,j] / ratio

D_p = simplify(D_p)
D_p


# #### Stiffness tensor

# In[79]:


C_isotropic_alpha_p = MutableDenseNDimArray.zeros(3, 3, 3, 3)
q=1+alpha3/R
for i in range(3):
    for j in range(3):        
        for k in range(3):
            for l in range(3):
                fact = 1
                if (i==0):
                    fact = fact*q
                if (j==0):
                    fact = fact*q
                if (k==0):
                    fact = fact*q
                if (l==0):
                    fact = fact*q
                C_isotropic_alpha_p[i,j,k,l] = simplify(C_isotropic_alpha[i,j,k,l]*fact)
            
C_isotropic_matrix_alpha_p = zeros(6)

for s in range(6):
    for t in range(6):
        i,j = getCIndecies(s)
        k,l = getCIndecies(t)
        C_isotropic_matrix_alpha_p[s,t] = C_isotropic_alpha_p[i,j,k,l]
        
C_isotropic_matrix_alpha_p


# In[80]:


C_orthotropic_alpha_p = MutableDenseNDimArray.zeros(3, 3, 3, 3)
q=1+alpha3/R
for i in range(3):
    for j in range(3):        
        for k in range(3):
            for l in range(3):
                fact = 1
                if (i==0):
                    fact = fact*q
                if (j==0):
                    fact = fact*q
                if (k==0):
                    fact = fact*q
                if (l==0):
                    fact = fact*q
                C_orthotropic_alpha_p[i,j,k,l] = simplify(C_orthotropic_alpha[i,j,k,l]*fact)
            
C_orthotropic_matrix_alpha_p = zeros(6)

for s in range(6):
    for t in range(6):
        i,j = getCIndecies(s)
        k,l = getCIndecies(t)
        C_orthotropic_matrix_alpha_p[s,t] = C_orthotropic_alpha_p[i,j,k,l]
        
C_orthotropic_matrix_alpha_p



# ## Square of segment 

# $A=\frac {\theta}{2} \left( R + h_2 \right)^2-\frac {\theta}{2} \left( R + h_1 \right)^2$

# In[82]:


theta, h1, h2=symbols('theta h_1 h_2')
square_geom=theta/2*(R+h2)**2-theta/2*(R+h1)**2
expand(simplify(square_geom))


# ${\displaystyle A=\int_{0}^{L}\int_{h_1}^{h_2} \left( 1+\frac{\alpha_3}{R} \right) d \alpha_1 d \alpha_3}, L=R \theta$

# In[83]:


square_int=integrate(integrate(1+alpha3/R, (alpha3, h1, h2)), (alpha1, 0, theta*R))
expand(simplify(square_int))


# ## Virtual work

# ### Isotropic material physical coordinates

# In[84]:

K=Symbol('K')
S = simplify(D_p.T*C_isotropic_matrix_alpha_p*D_p*(1+alpha3/R))
S = simplify(S.subs(R, 1/K))
S


# ## Mass matrix in physical coordinates

# In[88]:


rho=Symbol('rho')
B_h=zeros(3,12)
B_h[0,0]=1
B_h[1,4]=1
B_h[1,8]=1
M=simplify(rho*P.T*B_h.T*G_con*B_h*P)
M


# In[89]:


M_p = rho*B_h.T*B_h*(1+alpha3/R)
M_p = simplify(M_p.subs(R, 1/K))
M_p


# In[90]:


mass_matrix_func = lambdify([K, rho, alpha3], M_p, "numpy")
mass_matrix_func(100, 200, 400)



stiffness_matrix_func = lambdify([K, mu, la, alpha3], S, "numpy")
stiffness_matrix_func(100, 200, 300, 400)


# In[95]:


import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.solver as s
import fem.mesh as me
import plot


def generate_layers(thickness, layers_count, material):
    layer_top = thickness / 2
    layer_thickness = thickness / layers_count
    layers = set()
    for i in range(layers_count):
        layer = m.Layer(layer_top - layer_thickness, layer_top, material, i)
        layers.add(layer)
        layer_top -= layer_thickness
    return layers


def solve(width, curvature, thickness, N, M):
    layers_count = 1
    layers = generate_layers(thickness, layers_count, mat.IsotropicMaterial.steel())
    mesh = me.Mesh.generate(width, layers, N, M, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
    
    print(mesh.nodes_count())
#    geometry = g.CorrugatedCylindricalPlate(width, curvature, corrugation_amplitude, corrugation_frequency)
    geometry = g.CylindricalPlate(width, curvature)
#    geometry = g.Geometry()
    model = m.Model(geometry, layers, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
    return s.solve(model, mesh, stiffness_matrix, mass_matrix)

def stiffness_matrix(material, geometry, x1, x2, x3):
    return stiffness_matrix_func(1/geometry.curvature, material.mu(), material.lam(), x3)

def mass_matrix(material, geometry, x1, x2, x3):
    return mass_matrix_func(1/geometry.curvature, material.rho, x3)



# r=2
# width = r*2*3.14
# curvature = 1/r

width = 2
curvature = 0.8
thickness = 0.05

N = 100
M = 4

def to_cartesian(x1, x2, x3):
    x=alpha_x(1/curvature, width, x1, x3)
    z=alpha_z(1/curvature, width, x1, x3)
    return x,0,z

results = solve(width, curvature, thickness, N, M)
results_index = 0

#print(results[results_index].mesh.elements)
#plot.plot_mesh(results[results_index].mesh, width, thickness)
plot.plot_init_and_deformed_geometry_in_cartesian(results[results_index], 0, width, -thickness / 2, thickness / 2, 0, to_cartesian)
#plot.plot_init_geometry(results[results_index].geometry, 0, width, -thickness / 2, thickness / 2, 0)
# plot.plot_strain(results[results_index], 0, width, -thickness / 2, thickness / 2, 0)


to_print = 20
if (len(results) < to_print):
    to_print = len(results)

for i in range(to_print):
    print(results[i].freq)

