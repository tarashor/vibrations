# -*- coding: utf-8 -*-
import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.mesh as me
import numpy as np
#import plot

import utils

import fem.shells1D.secondorder.shellsolver as s1D2
import fem.shells1D.secondorder.result1D as r1D2
import fem.shells1D.secondorder.matrices1D as matrix1D2


def solve1D2(geometry, thickness, material, N):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
    mesh = me.Mesh.generate1D(geometry.width, layers, N, model.boundary_conditions)
    
    lam, vec = s1D2.solve(model, mesh, matrix1D2.stiffness_matrix, matrix1D2.mass_matrix)
    
    results = r1D2.Result.convert_to_results(lam, vec, mesh, geometry, thickness)
    
    return results

E = 40000
v = 0.25
rho = 2000

k = 40


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

results1D2 = solve1D2(geometry, thickness, material, N)
    
print("[1D2] w_min1 = {}".format(results1D2[0].freqHz()))
print("[1D2] w_min2 = {}".format(results1D2[1].freqHz()))
print("[1D2] w_min3 = {}".format(results1D2[2].freqHz()))

