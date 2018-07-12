import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.shells1D.firstorder.shellsolver as s
import fem.shells1D.firstorder.result1D as r

import fem.shells1D.firstorder.nonlinearshellsolver as s_nl
import fem.shells1D.firstorder.result1Dnonlinear as r_nl


import fem.mesh as me
import os
import platform
import matplotlib.pyplot as plt

import numpy as np

from fem.shells1D.firstorder.matrices1D import stiffness_matrix, mass_matrix, stiffness_matrix_nl_1, stiffness_matrix_nl_2


def solveLinear(geometry, thickness, material, N, bc):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, bc)
#    model = m.Model(geometry, layers, m.Model.FIXED_LEFT_RIGHT_EDGE)
    mesh = me.Mesh.generate1D(geometry.width, layers, N, model.boundary_conditions)
    
    lam, vec = s.solve(model, mesh, stiffness_matrix, mass_matrix)
    
    results_index = 0
    results = r.Result.convert_to_results(lam, vec, mesh, geometry)
    
    return results[results_index]


def solveNonlinear(geometry, thickness, material, N, u_max, bc):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, bc)
#    model = m.Model(geometry, layers, m.Model.FIXED_LEFT_RIGHT_EDGE)
    mesh = me.Mesh.generate1D(geometry.width, layers, N, model.boundary_conditions)
    
    lam_nl, res, U1, U2, U3, n = s_nl.solve_nl(model, mesh, stiffness_matrix, mass_matrix, stiffness_matrix_nl_1, stiffness_matrix_nl_2, u_max)
    
    return r_nl.ResultNL.convert_to_result(lam_nl, res, mesh, geometry), n



E = 40*(10**9)
#E = 40000
v = 0.3
rho = 8000

#kE3 = 100000000
kE3 = 1
#kG13 = 100000000
kG13 = 1

material = mat.OrthotropicMaterial.create_from_E_and_v_with_koef_E3(E,v,rho)

material.C[2,2] *= kE3
material.C[4,4] *= kG13

width = 1
curvature = 0
thickness = 0.01

corrugation_amplitude = 0
corrugation_frequency = 0

#bc = m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS

bc = m.Model.FIXED_LEFT_RIGHT_EDGE

geometry = g.General(width, curvature, corrugation_amplitude, corrugation_frequency)

N = 100

norm_koef = 0.2

result = solveLinear(geometry, thickness, material, N, bc)

x = []
y = []

yv = []

K = 3/4

if (bc == m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS):
    K = 3

for i in range(20):
    u_max = i*norm_koef*thickness
    resultNl, n = solveNonlinear(geometry, thickness, material, N, u_max, bc)
#    resultNl2 = solveNonlinear2(geometry, thickness, material, N, M, u_max)
    
#    print('w_max = {}, w_l = {}, w_nl = {}'.format(u_max,result.freqHz(), resultNl.freqHz()))
    
    print(n)
    
    d = i*norm_koef/n
    dy = resultNl.freqHz()/result.freqHz()
    dya = np.sqrt(1+0.75*K*d*d)
    x.append(d)
    y.append(dy)
    yv.append(dya)
    
    print('fem {} = {}'.format(d, dy))
    print('anal {} = {}'.format(d, dya))
    
    
tex_path = '/usr/local/texlive/2017/bin/x86_64-darwin'
if (platform.system() == 'Windows'):
    tex_path = "C:\Program Files\MiKTeX 2.9\miktex/bin/x64"

os.environ["PATH"] += os.pathsep + tex_path

plt.rc('text', usetex=True)
   
plt.rc('font', family='serif')

SMALL_SIZE = 24
MEDIUM_SIZE = 28
BIGGER_SIZE = 32

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

#plt.plot(x, y, 'o-', linewidth=2.0, markersize=8, markeredgewidth=2, markeredgecolor='r', markerfacecolor='None', label = "General 2D theory")
#plt.plot(x1D1, y1D1, 'x--', linewidth=2.0, markersize=8, markeredgewidth=2, markeredgecolor='b', markerfacecolor='None', label = "Mindlin-Reissner theory")
plt.plot(y, x, 'ro-', linewidth=2.0, markersize=8, markeredgewidth=2, markeredgecolor='r', markerfacecolor='None', label = "Shell second order theory")

plt.plot(yv, x, 'bv:', linewidth=2.0, markersize=8, markeredgewidth=2, markeredgecolor='r', markerfacecolor='None', label = "Volmir")
plt.xlabel(r"$\frac{\omega_{NL}}{\omega_{L}}$")
plt.ylabel(r"$\frac{w_{max}}{h}$")

plt.xlim(xmin=0)
plt.legend(loc='best')

#plt.title(r"Shear influence")
#plt.title(r"Convergence $\frac{E}{E_3} = 1$")
#plt.title("Збіжність")

plt.grid()
plt.show()
    
    

