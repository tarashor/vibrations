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
x = (R + ga * cos(gv * a1)) * cos(a1)
y = alpha2
z = (R + ga * cos(gv * a1)) * sin(a1)
r = x*N.i + y*N.j + z*N.k

r1 = trigsimp(r.diff(alpha1))
r2 = trigsimp(r.diff(alpha2))


#r1m=trigsimp(simplify(collect(expand(r1.dot(r1)), gv)))
#printing.latex.print_latex(r1m)

r1m=sympify("((R + g_a*cos(g_v*(L + pi*R - 2*alpha_1)/(2*R)))**2 + g_a**2*g_v**2*sin(g_v*(L + pi*R - 2*alpha_1)/(2*R))**2)/R**2")

display(r1m)

k1 = r1m**(S(1)/S(2))
k2 = trigsimp(r2.magnitude())
r1 = r1/k1
r2 = r2/k2

display(r1)
##
##
n = r1.cross(r2)
##
n_len = trigsimp(n.dot(n))
##
display(n_len)


