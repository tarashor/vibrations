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
import plot

import utils


def solveLinear2D(geometry, layers, N, M, bc):
    model = m.Model(geometry, layers, bc)
    mesh = me.Mesh.generate2D(geometry.width, layers, N, M, model.boundary_conditions)
    
    lam, vec = s2D.solve(model, mesh, mat2D.stiffness_matrix, mat2D.mass_matrix)
    
    results_index = 0
    results = r2D.Result.convert_to_results(lam, vec, mesh, geometry)
    
#    for i in range(2):
#        plot.plot_init_and_deformed_geometry_in_cartesian(results[i], 0, width, -thickness / 2, thickness / 2, 0, geometry.to_cartesian_coordinates)
    
    return results[results_index]


def solveNonlinear2D(geometry, layers, N, M, u_max, bc):
    model = m.Model(geometry, layers, bc)
    mesh = me.Mesh.generate2D(geometry.width, layers, N, M, model.boundary_conditions)
    
    lam_nl, res, U1, U2, U3, n = s2Dnl.solve_nl(model, mesh, mat2D.stiffness_matrix, mat2D.mass_matrix, mat2D.stiffness_matrix_nl_1, mat2D.stiffness_matrix_nl_2, u_max)
    
    return r2Dnl.ResultNL.convert_to_result(lam_nl, res, U1, U2, U3, mesh, geometry), n



width = 1
thickness = 0.1

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
M = 2

bc = m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS
#bc = m.Model.FIXED_LEFT_RIGHT_EDGE


curvatures = [0.5, 0.8, 1, 1.5, 2]
#curvatures = [0.05]

y_per_hr = {}
x_per_hr = {}

for curvature in curvatures:

    geometry = g.General(width, curvature, corrugation_amplitude, corrugation_frequency)
    
    norm_koef = 0.5
    
    result2D = solveLinear2D(geometry, layers, N, M, bc)
    
    
    print('========================')
    print(curvature)
    print(result2D.freqHz())
    print('========================')
    
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
        
        
    y_per_hr[curvature] = y2D
    x_per_hr[curvature] = x
        

folder = "./results/curved/"

utils.save_results(folder+"y_per_hr", y_per_hr)
utils.save_results(folder+"x_per_hr", x_per_hr)

    
    

