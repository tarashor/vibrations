import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.mesh as me


import fem.shells1D.secondorder.shellsolver as s1D2O
import fem.shells1D.secondorder.result1D as r1D2O
import fem.shells1D.secondorder.matrices1D as mat1D2O
import fem.shells1D.secondorder.nonlinearshellsolver as s1D2Onl
import fem.shells1D.secondorder.result1Dnonlinear as r1D2Onl

import os
import platform
import matplotlib.pyplot as plt

import utils


def solveLinear1D2O(geometry, layers, N, bc):
    model = m.Model(geometry, layers, bc)
#    model = m.Model(geometry, layers, m.Model.FIXED_LEFT_RIGHT_EDGE)
    mesh = me.Mesh.generate1D(geometry.width, layers, N, model.boundary_conditions)
    
    lam, vec = s1D2O.solve(model, mesh, mat1D2O.stiffness_matrix, mat1D2O.mass_matrix)
    
    results_index = 0
    results = r1D2O.Result.convert_to_results(lam, vec, mesh, geometry, thickness)
    
    return results[results_index]


def solveNonlinear1D2O(geometry, layers, N, u_max, bc):

    model = m.Model(geometry, layers, bc)
#    model = m.Model(geometry, layers, m.Model.FIXED_LEFT_RIGHT_EDGE)
    mesh = me.Mesh.generate1D(geometry.width, layers, N, model.boundary_conditions)
    
    lam_nl, res, U1, U2, U3, n = s1D2Onl.solve_nl(model, mesh, mat1D2O.stiffness_matrix, mat1D2O.mass_matrix, mat1D2O.stiffness_matrix_nl_1, mat1D2O.stiffness_matrix_nl_2, u_max)
    
    return r1D2Onl.ResultNL.convert_to_result(lam_nl, res, mesh, geometry, thickness), n



width = 1
thickness = 0.01

corrugation_amplitude = 0
corrugation_frequency = 0

E = 40*(10**9)
v = 0.25
rho = 8000

#kE3 = 100000000
kE1 = 40
#kG13 = 100000000
kG13 = 1

material = mat.OrthotropicMaterial.create_from_E_and_v_with_koef_E1(E,v,rho, kE1)

material.C[4,4] *= kG13


layers = m.Layer.generate_layers(thickness, [material])

N = 50
M = 4

bc = m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS
#bc = m.Model.FIXED_LEFT_RIGHT_EDGE


curvatures = [0, 0.5, 0.8, 1, 1.5, 2]

y_per_hr = {}
x_per_hr = {}

for curvature in curvatures:

    geometry = g.General(width, curvature, corrugation_amplitude, corrugation_frequency)
    
    norm_koef = 0.5
    
    result2D = solveLinear1D2O(geometry, layers, N, bc)
    
    print('========================')
    print(curvature)
    print(result2D.freqHz())
    print('========================')
    
    x = []
    y2D = []
    
    
    for i in range(4):
        u_max = i*norm_koef*thickness
        result2Dnl, n = solveNonlinear1D2O(geometry, layers, N, u_max, bc)
        
        
        d = i*norm_koef
        dy2D = result2Dnl.freqHz()/result2D.freqHz()
        
        x.append(d)
        y2D.append(dy2D)
        
        print('2D {} = {}'.format(d, dy2D))
        
        
    y_per_hr[curvature] = y2D
    x_per_hr[curvature] = x
        

folder = "./results/curved/"

utils.save_results(folder+"y_per_hr", y_per_hr)
utils.save_results(folder+"x_per_hr", x_per_hr)

    
    

