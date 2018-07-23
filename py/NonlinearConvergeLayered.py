import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.mesh as me


import fem.general2D.solverlinear as s2D
import fem.general2D.result2D as r2D
import fem.general2D.matrices2D as mat2D
import fem.general2D.solver_nonlinear as s2Dnl
import fem.general2D.result2Dnonlinear as r2Dnl

import os
import platform
import matplotlib.pyplot as plt

import numpy as np

import plot as p

def solveLinear2D(geometry, layers, N, M, bc):
    model = m.Model(geometry, layers, bc)
    mesh = me.Mesh.generate2D(geometry.width, layers, N, M, model.boundary_conditions)
    
    lam, vec = s2D.solve(model, mesh, mat2D.stiffness_matrix, mat2D.mass_matrix)
    
    results_index = 0
    results = r2D.Result.convert_to_results(lam, vec, mesh, geometry)
    
    return results[results_index]


def solveNonlinear2D(geometry, layers, N, M, u_max, bc):
    model = m.Model(geometry, layers, bc)
    mesh = me.Mesh.generate2D(geometry.width, layers, N, M, model.boundary_conditions)
#    p.plot_mesh(mesh, 1, 0.1)
    
    lam_nl, res, U1, U2, U3, n = s2Dnl.solve_nl(model, mesh, mat2D.stiffness_matrix, mat2D.mass_matrix, mat2D.stiffness_matrix_nl_1, mat2D.stiffness_matrix_nl_2, u_max)
    
    return r2Dnl.ResultNL.convert_to_result(lam_nl, res, U1, U2, U3, mesh, geometry), n





width = 1
curvature = 0
thickness = 0.1

corrugation_amplitude = 0
corrugation_frequency = 0

geometry = g.General(width, curvature, corrugation_amplitude, corrugation_frequency)

N = 50
M = 4

bc = m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS

#bc = m.Model.FIXED_LEFT_RIGHT_EDGE

material1 = mat.IsotropicMaterial.steel()
material2 = mat.IsotropicMaterial.rubber()

rubber_coefs = [0, 0.4, 0.6, 0.8, 0.9, 0.95, 1]

#rubber_coefs = [0.2, 0.4, 0.6]

y_per_hr = {}
x_per_hr = {}

for rubber_h_coef in rubber_coefs:

    layers = set()
    if (rubber_h_coef != 0 and rubber_h_coef != 1):
        rh = rubber_h_coef*thickness
        sh = (1-rubber_h_coef)*thickness/2
        
        h = sh+rh+sh
        top = h/2
        layer1 = m.Layer(top - sh, top, material1, 0)
        layer2 = m.Layer(top - sh-rh, top-sh, material2, 1)
        layer3 = m.Layer(top - sh-rh - sh, top - sh-rh, material1, 2)
        layers.add(layer1)
        layers.add(layer2)
        layers.add(layer3)
    elif (rubber_h_coef == 0):
        sh = thickness
        
        top = thickness/2
        layer1 = m.Layer(top - sh, top, material1, 0)
        
        layers.add(layer1)
        
    elif (rubber_h_coef == 1):
        rh = thickness
        
        top = thickness/2
        layer1 = m.Layer(top - rh, top, material2, 0)
        
        layers.add(layer1)
        
        
    
    
    norm_koef = 0.5
    
    result2D = solveLinear2D(geometry, layers, N, M, bc)
    
    
    x = []
    
    y2D = []
    
    
    for i in range(7):
        u_max = i*norm_koef*thickness
        result2Dnl, n = solveNonlinear2D(geometry, layers, N, M, u_max, bc)
        
        
        d = i*norm_koef
        dy2D = result2Dnl.freqHz()/result2D.freqHz()
        
        
        x.append(d)
        y2D.append(dy2D)
        
        print('2D {} = {}'.format(d, dy2D))
        
        
    y_per_hr[rubber_h_coef] = y2D
    x_per_hr[rubber_h_coef] = x
        
        
    
plt.figure()
    
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

#plt.plot(y2D, x, 'ro-', linewidth=2.0, markersize=7, markeredgewidth=2, markeredgecolor='r', markerfacecolor='r', label = "Total nonlinearity")
#plt.plot(y1D2O, x, 'gs--', linewidth=2.0, markersize=7, markeredgewidth=2, markeredgecolor='g', markerfacecolor='g', label = "Second order")
#plt.plot(y1D1O, x, 'bx-.', linewidth=2.0, markersize=7, markeredgewidth=2, markeredgecolor='b', markerfacecolor='b', label = "Mindlin-Reisner")

#plt.plot(ya, x, 'mv:', linewidth=2.0, markersize=7, markeredgewidth=2, markeredgecolor='m', markerfacecolor='m', label = "Analytical")

for hr in rubber_coefs:
    print('===========')
    print(hr)
    print(y_per_hr[hr])
    print(x_per_hr[hr])
    
    plt.plot(y_per_hr[hr], x_per_hr[hr], label = r"$h_r = {}$".format(hr))

plt.xlabel(r"$\frac{\omega_{NL}}{\omega_{L}}$")
plt.ylabel(r"$\frac{w_{max}}{h}$")

plt.xlim(xmin=0)
plt.legend(loc='best')

plt.grid()
plt.show()
    
    

