# -*- coding: utf-8 -*-
import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.mesh as me
import numpy as np
#import plot

import utils

import fem.general2D.solverlinear as s
import fem.general2D.result2D as r

import plot

from fem.general2D.matrices2D import stiffness_matrix, mass_matrix, stiffness_matrix_nl_1, stiffness_matrix_nl_2




def solve(geometry, thickness, material, N, M):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
#    model = m.Model(geometry, layers, m.Model.FIXED_LEFT_RIGHT_EDGE)
    mesh = me.Mesh.generate2D(geometry.width, layers, N, M, model.boundary_conditions)
    
    lam, vec = s.solve(model, mesh, stiffness_matrix, mass_matrix)
    
    return lam, vec, mesh, geometry


#E = 10*(10**9)
E=1
v = 0.25
rho = 2000

k = 3


Es = [k*E, E, E]
vs = np.zeros((3,3))

for i in range(3):
    for j in range(3):
        if (i != j):
            vs[i,j] = v
            
vs[2,0] = (v/k)

G = 0.5*E
Gs = [G,G,G]

material = mat.OrthotropicMaterial.create_from_E_and_v(Es,vs, Gs,rho)

#values_N_width = [10, 25, 50, 75, 100, 150, 200, 300]

width = 1
curvature = 0
thickness = 0.01


corrugation_amplitude = 0
corrugation_frequency = 0

geometry = g.General(width, curvature, corrugation_amplitude, corrugation_frequency)

N = 50
M = 4

lam, vec, mesh, geometry = solve(geometry, thickness, material, N, M)
results = r.Result.convert_to_results(lam, vec, mesh, geometry)

print("w_min1 = {}".format(results[0].freqHz()))
print("w_min2 = {}".format(results[1].freqHz()))
print("w_min3 = {}".format(results[2].freqHz()))

