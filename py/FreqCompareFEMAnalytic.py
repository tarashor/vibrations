import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.mesh as me

import fem.general2D.solverlinear as s
import fem.general2D.result2D as r

import fem.general2D.solver_nonlinear as s_nl
import fem.general2D.result2Dnonlinear as r_nl

import plot

import numpy as np

from fem.general2D.matrices2D import stiffness_matrix, mass_matrix, stiffness_matrix_nl_1, stiffness_matrix_nl_2


def solveLinear(geometry, thickness, material, N, M, u_max):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
#    model = m.Model(geometry, layers, m.Model.FIXED_LEFT_RIGHT_EDGE)
    mesh = me.Mesh.generate2D(geometry.width, layers, N, M, model.boundary_conditions)
    
    lam, vec = s.solve(model, mesh, stiffness_matrix, mass_matrix, u_max)
    
    results_index = 0
    results = r.Result.convert_to_results(lam, vec, mesh, geometry)
    
    return results[results_index]


def wAnalyticalLin(geometry, thickness, material, N, M, u_max):
    
    ko = 14/15
    G = material.mu()
    LAM = ko * thickness * G
    
    c2 = np.sqrt(LAM/(rho*thickness))
    
    lam = np.pi / geometry.width
    
    alpha = (1+v)*v*v/(1-v-2*v*v)
    
    B = E*thickness/(1-v*v)*(1+alpha)
    
    D = (thickness*thickness/12)*B
    
    k12 = LAM / D 
    
    w = c2*lam*lam/(np.sqrt(k12+lam*lam))
    
    return w




E = 40000
#E = 40000
v = 0.3
rho = 2000

material = mat.IsotropicMaterial(E,v,rho)

width = 1
curvature = 0
thickness = 0.1

corrugation_amplitude = 0
corrugation_frequency = 0

geometry = g.General(width, curvature, corrugation_amplitude, corrugation_frequency)

N = 100
M = 4

norm_koef = 2
u_max = norm_koef*thickness

result = solveLinear(geometry, thickness, material, N, M, u_max)

wa = wAnalyticalLin(geometry, thickness, material, N, M, u_max)

print('FEM = {}'.format(result.freqHz()))

print('Analyt = {}'.format(wa))
