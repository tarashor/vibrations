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

Ralpha = r+alpha3*n

R1=Ralpha.diff(alpha1)
R2=Ralpha.diff(alpha2)
R3=Ralpha.diff(alpha3)

display(R1)
display(R3)

