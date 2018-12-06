import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.mesh as me


import fem.shells1D.secondorder.shellsolver as s1D2O
import fem.shells1D.secondorder.result1D as r1D2O
import fem.shells1D.secondorder.matrices1D as mat1D2O

import fem.general2D.solverlinear as s2D
import fem.general2D.result2D as r2D
import fem.general2D.matrices2D as mat2D

import plot

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter


def solveLinear1D2O(geometry, layers, N, bc):
    model = m.Model(geometry, layers, bc)
#    model = m.Model(geometry, layers, m.Model.FIXED_LEFT_RIGHT_EDGE)
    mesh = me.Mesh.generate1D(geometry.width, layers, N, model.boundary_conditions)
    
    lam, vec = s1D2O.solve(model, mesh, mat1D2O.stiffness_matrix, mat1D2O.mass_matrix)
    
    results_index = 0
    results = r1D2O.Result.convert_to_results(lam, vec, mesh, geometry, thickness)
    
    return results[results_index]

def solveLinear2D(geometry, layers, N, M, bc):
    model = m.Model(geometry, layers, bc)
    mesh = me.Mesh.generate2D(geometry.width, layers, N, M, model.boundary_conditions)
    
    lam, vec = s2D.solve(model, mesh, mat2D.stiffness_matrix, mat2D.mass_matrix)
    
    results_index = 0
    results = r2D.Result.convert_to_results(lam, vec, mesh, geometry)
    
    return results[results_index]

width = 2
curvature = 1.25
thickness = 0.05

corrugation_amplitudes = [0, 0.015, 0.03, 0.06, 0.1, 0.15, 0.2, 0.25, 0.3]
corrugation_frequencies = [2, 6, 15, 20, 50, 75, 100]

E = 40*(10**9)
v = 0.25
rho = 8000

#kE3 = 100000000
kE1 = 40
#kG13 = 100000000
kG13 = 1

#material = mat.OrthotropicMaterial.create_from_E_and_v_with_koef_E1(E,v,rho, kE1)

#material.C[4,4] *= kG13

material = mat.IsotropicMaterial.steel()

N = 50

#bc = m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS
bc = m.Model.FIXED_LEFT_RIGHT_EDGE


v = np.zeros((len(corrugation_amplitudes), len(corrugation_frequencies)))

layers = m.Layer.generate_layers(thickness, [material])

i= 0

for corrugation_amplitude in corrugation_amplitudes:
    j = 0
    for corrugation_frequency in corrugation_frequencies:
        
        geometry = g.General(width, curvature, corrugation_amplitude, corrugation_frequency)
        
        print("g_a =  {}".format(corrugation_amplitude))
        print("g_v =  {}".format(corrugation_frequency))
        
        result2D = solveLinear1D2O(geometry, layers, N, bc)
#        result2D = solveLinear2D(geometry, layers, N, 2, bc)
        
        v[i, j] = result2D.freqHz()
        
        j+=1
    
    i+=1
        
plot.init()
fig=plt.figure()

(X1, X2) = np.meshgrid(corrugation_amplitudes, corrugation_frequencies)
surf = plt.contourf(X1, X2, v.T, 100)
#surf = plt.contourf(X1, X2, v.T)
plt.ylabel(r"$g_v$")
plt.xlabel(r"$g_A, m$")
#plt.clabel(surf, inline=1, fontsize=16)
plt.colorbar(surf)
plt.grid()
plt.show()

#ax = fig.gca(projection='3d')
#surf = ax.plot_surface(X1, X2, v.T, cmap=cm.rainbow)
#ax.yaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#plt.show()
        
    
        
    
    

