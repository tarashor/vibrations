import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.mesh as me

import fem.shells1D.firstorder.shellsolver as s1D1O
import fem.shells1D.firstorder.result1D as r1D1O
import fem.shells1D.firstorder.matrices1D as mat1D1O
import fem.shells1D.firstorder.nonlinearshellsolver as s1D1Onl
import fem.shells1D.firstorder.result1Dnonlinear as r1D1Onl

import fem.shells1D.secondorder.shellsolver as s1D2O
import fem.shells1D.secondorder.result1D as r1D2O
import fem.shells1D.secondorder.matrices1D as mat1D2O
import fem.shells1D.secondorder.nonlinearshellsolver as s1D2Onl
import fem.shells1D.secondorder.result1Dnonlinear as r1D2Onl

import fem.general2D.solverlinear as s2D
import fem.general2D.result2D as r2D
import fem.general2D.matrices2D as mat2D
import fem.general2D.solver_nonlinear as s2Dnl
import fem.general2D.result2Dnonlinear as r2Dnl

import os
import platform
import matplotlib.pyplot as plt

import numpy as np

def solveLinear2D(geometry, thickness, material, N, M, bc):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, bc)
    mesh = me.Mesh.generate2D(geometry.width, layers, N, M, model.boundary_conditions)
    
    lam, vec = s2D.solve(model, mesh, mat2D.stiffness_matrix, mat2D.mass_matrix)
    
    results_index = 0
    results = r2D.Result.convert_to_results(lam, vec, mesh, geometry)
    
    return results[results_index]


def solveNonlinear2D(geometry, thickness, material, N, M, u_max, bc):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, bc)
    mesh = me.Mesh.generate2D(geometry.width, layers, N, M, model.boundary_conditions)
    
    lam_nl, res, U1, U2, U3, n = s2Dnl.solve_nl(model, mesh, mat2D.stiffness_matrix, mat2D.mass_matrix, mat2D.stiffness_matrix_nl_1, mat2D.stiffness_matrix_nl_2, u_max)
    
    return r2Dnl.ResultNL.convert_to_result(lam_nl, res, U1, U2, U3, mesh, geometry), n

def solveLinear1D1O(geometry, thickness, material, N, bc):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, bc)
    mesh = me.Mesh.generate1D(geometry.width, layers, N, model.boundary_conditions)
    
    lam, vec = s1D1O.solve(model, mesh, mat1D1O.stiffness_matrix, mat1D1O.mass_matrix)
    
    results_index = 0
    results = r1D1O.Result.convert_to_results(lam, vec, mesh, geometry)
    
    return results[results_index]


def solveNonlinear1D1O(geometry, thickness, material, N, u_max, bc):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, bc)
#    model = m.Model(geometry, layers, m.Model.FIXED_LEFT_RIGHT_EDGE)
    mesh = me.Mesh.generate1D(geometry.width, layers, N, model.boundary_conditions)
    
    lam_nl, res, U1, U2, U3, n = s1D1Onl.solve_nl(model, mesh, mat1D1O.stiffness_matrix, mat1D1O.mass_matrix, mat1D1O.stiffness_matrix_nl_1, mat1D1O.stiffness_matrix_nl_2, u_max)
    
    return r1D1Onl.ResultNL.convert_to_result(lam_nl, res, mesh, geometry), n

def solveLinear1D2O(geometry, thickness, material, N, bc):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, bc)
#    model = m.Model(geometry, layers, m.Model.FIXED_LEFT_RIGHT_EDGE)
    mesh = me.Mesh.generate1D(geometry.width, layers, N, model.boundary_conditions)
    
    lam, vec = s1D2O.solve(model, mesh, mat1D2O.stiffness_matrix, mat1D2O.mass_matrix)
    
    results_index = 0
    results = r1D2O.Result.convert_to_results(lam, vec, mesh, geometry, thickness)
    
    return results[results_index]


def solveNonlinear1D2O(geometry, thickness, material, N, u_max, bc):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, bc)
#    model = m.Model(geometry, layers, m.Model.FIXED_LEFT_RIGHT_EDGE)
    mesh = me.Mesh.generate1D(geometry.width, layers, N, model.boundary_conditions)
    
    lam_nl, res, U1, U2, U3, n = s1D2Onl.solve_nl(model, mesh, mat1D2O.stiffness_matrix, mat1D2O.mass_matrix, mat1D2O.stiffness_matrix_nl_1, mat1D2O.stiffness_matrix_nl_2, u_max)
    
    return r1D2Onl.ResultNL.convert_to_result(lam_nl, res, mesh, geometry, thickness), n


def getK(geometry, thickness, material, bc):
    K = 3/4

    if (bc == m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS):
        K = 3
        
    l = thickness/geometry.width
    
    K += np.pi*np.pi*l*l*material.C[0,0]/material.C[4,4]/4
        
    return K

E = 40*(10**9)
#E = 40000
v = 0.25
rho = 8000

#kE3 = 100000000
kE1 = 40
#kG13 = 100000000
kG13 = 1

material = mat.OrthotropicMaterial.create_from_E_and_v_with_koef_E1(E,v,rho, kE1)

material.C[4,4] *= kG13

width = 1
curvature = 0
thickness = 0.01

corrugation_amplitude = 0
corrugation_frequency = 0

geometry = g.General(width, curvature, corrugation_amplitude, corrugation_frequency)

N = 100
M = 4

bc = m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS

#bc = m.Model.FIXED_LEFT_RIGHT_EDGE


norm_koef = 0.2

result2D = solveLinear2D(geometry, thickness, material, N, M, bc)
result1D2O = solveLinear1D2O(geometry, thickness, material, N, bc)
result1D1O = solveLinear1D1O(geometry, thickness, material, N, bc)

x = []

y2D = []
y1D2O = []
y1D1O = []
ya = []

K = getK(geometry, thickness, material, bc)

print(K)

for i in range(7):
    u_max = i*norm_koef*thickness
    result2Dnl, n = solveNonlinear2D(geometry, thickness, material, N, M, u_max, bc)
    result1D1Onl, n = solveNonlinear1D1O(geometry, thickness, material, N, u_max, bc)
    result1D2Onl, n = solveNonlinear1D2O(geometry, thickness, material, N, u_max, bc)
    
    d = i*norm_koef
    dy2D = result2Dnl.freqHz()/result2D.freqHz()
    dy1D2O = result1D2Onl.freqHz()/result1D2O.freqHz()
    dy1D1O = result1D1Onl.freqHz()/result1D1O.freqHz()
    
    dya = np.sqrt(1+0.75*K*d*d)
    
    x.append(d)
    y2D.append(dy2D)
    y1D2O.append(dy1D2O)
    y1D1O.append(dy1D1O)
    ya.append(dya)
    
    print('2D {} = {}'.format(d, dy2D))
    print('1D2O {} = {}'.format(d, dy1D2O))
    print('1D1O {} = {}'.format(d, dy1D1O))
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

plt.plot(y2D, x, 'ro-', linewidth=2.0, markersize=7, markeredgewidth=2, markeredgecolor='r', markerfacecolor='r', label = "Total nonlinearity")
plt.plot(y1D2O, x, 'gs--', linewidth=2.0, markersize=7, markeredgewidth=2, markeredgecolor='g', markerfacecolor='g', label = "Second order")
plt.plot(y1D1O, x, 'bx-.', linewidth=2.0, markersize=7, markeredgewidth=2, markeredgecolor='b', markerfacecolor='b', label = "Mindlin-Reisner")

plt.plot(ya, x, 'mv:', linewidth=2.0, markersize=7, markeredgewidth=2, markeredgecolor='m', markerfacecolor='m', label = "Analytical")
plt.xlabel(r"$\frac{\omega_{NL}}{\omega_{L}}$")
plt.ylabel(r"$\frac{w_{max}}{h}$")

plt.xlim(xmin=0)
plt.legend(loc='best')

#plt.title(r"Shear influence")
#plt.title(r"Convergence $\frac{E}{E_3} = 1$")
#plt.title("Збіжність")

plt.grid()
plt.show()
    
    

