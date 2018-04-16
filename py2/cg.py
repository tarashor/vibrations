#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 25 21:15:42 2018

@author: tarasgoriachko
"""

#import os
#
#os.environ["PATH"] += os.pathsep + '/usr/local/texlive/2017/bin/x86_64-darwin'

from IPython.display import display


from sympy import *
from sympy.vector import CoordSys3D
N = CoordSys3D('N')
x1, x2, x3 = symbols("x_1 x_2 x_3")
alpha1, alpha2, alpha3 = symbols("alpha_1 alpha_2 alpha3")
R, L, ga, gv = symbols("R L g_a g_v")
init_printing()

a1 = pi / 2 + (L / 2 - alpha1)/R

a2 = 2 * pi * alpha1 / L

x = (R + ga * cos(gv * a2)) * cos(a1)
y = alpha2
z = (R + ga * cos(gv * a2)) * sin(a1)
r = x*N.i + y*N.j + z*N.k
#r1 = trigsimp(r.diff(alpha1))

#r1 = ((-ga*gv*sin((L/2 - alpha_1)/R)*sin(g_v*(L + pi*R - 2*alpha_1)/(2*R)) + (R + g_a*cos(g_v*(L + pi*R - 2*alpha_1)/(2*R)))*cos((L/2 - alpha_1)/R))/sqrt(g_a**2*g_v**2*sin(g_v*(L + pi*R - 2*alpha_1)/(2*R))**2 + (R + g_a*cos(g_v*(L + pi*R - 2*alpha_1)/(2*R)))**2))*N.i + ((g_a*g_v*sin(g_v*(L + pi*R - 2*alpha_1)/(2*R))*cos((L/2 - alpha_1)/R) + (R + g_a*cos(g_v*(L + pi*R - 2*alpha_1)/(2*R)))*sin((L/2 - alpha_1)/R))/sqrt(g_a**2*g_v**2*sin(g_v*(L + pi*R - 2*alpha_1)/(2*R))**2 + (R + g_a*cos(g_v*(L + pi*R - 2*alpha_1)/(2*R)))**2))*N.k

r1 = (ga*gv*cos(a1)*sin(gv*a1) + (R + ga*cos(gv*a1))*sin(a1))*N.i + (ga*gv*sin(gv*a1)*sin(a1) - (R + ga*cos(gv*a1))*cos(a1))*N.k
r2 = trigsimp(r.diff(alpha2))

#r1m=r1.dot(r1)
r1m=(R + ga*cos(gv*a1))**2+(ga*gv*sin(gv*a1))**2
#r1m=sympify("((R + g_a*cos(g_v*(L + pi*R - 2*alpha_1)/(2*R)))**2 + g_a**2*g_v**2*sin(g_v*(L + pi*R - 2*alpha_1)/(2*R))**2)")

display(r1m)

k1 = sqrt(r1m)
r1 = simplify(r1/k1)


display(r1)

#dr1=r1.diff(alpha1)

dr1=((ga*gv**2*(-1)*cos(a1)*cos(gv*a1)/R + 2*ga*gv*sin(gv*a1)*sin(a1)/R + (R + ga*cos(gv*a1))*(-1)*cos(a1)/R)/sqrt(ga**2*gv**2*sin(gv*a1)**2 + (R + ga*cos(gv*a1))**2) + (-ga*gv*(-1)*cos(a1)*sin(gv*a1) + (R + ga*cos(gv*a1))*sin(a1))*(ga**2*gv**3*sin(gv*a1)*cos(gv*a1)/R - ga*gv*(R + ga*cos(gv*a1))*sin(gv*a1)/R)/(ga**2*gv**2*sin(gv*a1)**2 + (R + ga*cos(gv*a1))**2)**(3/2))*N.i + ((-ga*gv**2*sin(a1)*cos(gv*a1)/R + 2*ga*gv*(-1)*cos(a1)*sin(gv*a1)/R - (R + ga*cos(gv*a1))*sin(a1)/R)/sqrt(ga**2*gv**2*sin(gv*a1)**2 + (R + ga*cos(gv*a1))**2) + (ga*gv*sin(gv*a1)*sin(a1) + (R + ga*cos(gv*a1))*(-1)*cos(a1))*(ga**2*gv**3*sin(gv*a1)*cos(gv*a1)/R - ga*gv*(R + ga*cos(gv*a1))*sin(gv*a1)/R)/(ga**2*gv**2*sin(gv*a1)**2 + (R + ga*cos(gv*a1))**2)**(3/2))*N.k


display(dr1)


#r2 = r2/k2
#
#display(r1)
###
###
n = r1.cross(r2)
###
n_len = trigsimp(n.dot(n))
n=n/n_len
###
display(n)

dn1=n.diff(alpha1)

display(dn1)

Ralpha = r+alpha3*n
#
R1=dr1+alpha3*dn1
R2=Ralpha.diff(alpha2)
R3=n
#
display(R1)
display(R3)

#%%

eps=trigsimp(R1.dot(R2.cross(R3)))
R_1=(R2.cross(R3)/eps)
R_2=(R3.cross(R1)/eps)
R_3=(R1.cross(R2)/eps)

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