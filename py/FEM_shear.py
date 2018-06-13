# -*- coding: utf-8 -*-
import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.mesh as me
#import plot

import utils
#
#import fem.general2D.solverlinear as s2D
#import fem.general2D.result2D as r2D
#import fem.general2D.matrices2D as matrix2D
#
#
#def solve2D(geometry, thickness, material, N, M):
#    layers = m.Layer.generate_layers(thickness, [material])
#    model = m.Model(geometry, layers, m.Model.FIXED_LEFT_RIGHT_EDGE)
#    mesh = me.Mesh.generate2D(geometry.width, layers, N, M, model.boundary_conditions)
#    
#    lam, vec = s2D.solve(model, mesh, matrix2D.stiffness_matrix, matrix2D.mass_matrix)
#    
#    results = r2D.Result.convert_to_results(lam, vec, mesh, geometry)
#    
#    return results
#
#
#import fem.shells1D.firstorder.shellsolver as s1D1
#import fem.shells1D.firstorder.result1D as r1D1
#import fem.shells1D.firstorder.matrices1D as matrix1D1
#
#
#def solve1D1(geometry, thickness, material, N, M):
#    layers = m.Layer.generate_layers(thickness, [material])
#    model = m.Model(geometry, layers, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
#    mesh = me.Mesh.generate1D(geometry.width, layers, N, model.boundary_conditions)
#    
#    lam, vec = s1D1.solve(model, mesh, matrix1D1.stiffness_matrix, matrix1D1.mass_matrix)
#    
#    results = r1D1.Result.convert_to_results(lam, vec, mesh, geometry)
#    
#    return results

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

E = 40000000000
v = 0.3
rho = 2000

#values_N_width = [10, 25, 50, 75, 100, 150, 200, 300]

kShear = range(1,110, 3)

width = 1
curvature = 0
thickness = 0.1

N = 150

corrugation_amplitude = 0
corrugation_frequency = 0

geometry = g.General(width, curvature, corrugation_amplitude, corrugation_frequency)

#results2D_all_n = {}

#results1D1_all_n = {}

results1D2_all_n = {}

for k in kShear:
    print("===== k = {} =======".format(k))
    
    material = mat.IsotropicMaterial(E,v,rho)
    material.C[4,4] /= k
    
    results1D2 = solve1D2(geometry, thickness, material, N)
    
    print("[1D2] w_min = {}".format(results1D2[0].freqHz()))
    
    results1D2_all_n[k]=results1D2
    
    
folder = "./results/shear/"
#folder = ""
#utils.save_results(folder+"rNs2D", results2D_all_n)
#utils.save_results(folder+"rNs1D1", results1D1_all_n)
utils.save_results(folder+"rNs1D2", results1D2_all_n)