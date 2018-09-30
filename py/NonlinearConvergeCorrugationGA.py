import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.mesh as me


import fem.general2D.solverlinear as s2D
import fem.general2D.result2D as r2D
import fem.general2D.matrices2D as mat2D
import fem.general2D.solver_nonlinear as s2Dnl
import fem.general2D.result2Dnonlinear as r2Dnl

import fem.shells1D.secondorder.shellsolver as s1D2O
import fem.shells1D.secondorder.result1D as r1D2O
import fem.shells1D.secondorder.matrices1D as mat1D2O
import fem.shells1D.secondorder.nonlinearshellsolver as s1D2Onl
import fem.shells1D.secondorder.result1Dnonlinear as r1D2Onl

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


def solveLinear1D2O(geometry, layers, N, bc):
    model = m.Model(geometry, layers, bc)
    mesh = me.Mesh.generate1D(geometry.width, layers, N, model.boundary_conditions)
    
    lam, vec = s1D2O.solve(model, mesh, mat1D2O.stiffness_matrix, mat1D2O.mass_matrix)
    
    results_index = 0
    results = r1D2O.Result.convert_to_results(lam, vec, mesh, geometry, thickness)
    
#    for i in range(6):
#        plot.plot_init_and_deformed_geometry_in_cartesian(results[i], 0, width, -thickness / 2, thickness / 2, 0, geometry.to_cartesian_coordinates)
    
    return results[results_index]

def solveNonlinear1D2O(geometry, layers, N, u_max, bc):
    model = m.Model(geometry, layers, bc)
    mesh = me.Mesh.generate1D(geometry.width, layers, N, model.boundary_conditions)
    
    lam_nl, res, U1, U2, U3, n = s1D2Onl.solve_nl(model, mesh, mat1D2O.stiffness_matrix, mat1D2O.mass_matrix, mat1D2O.stiffness_matrix_nl_1, mat1D2O.stiffness_matrix_nl_2, u_max)
    
    return r1D2Onl.ResultNL.convert_to_result(lam_nl, res, mesh, geometry, thickness), n



width = 2
curvature = 0.8
thickness = 0.05
corrugation_frequency = 20

corrugation_amplitude = 0.03

material = mat.IsotropicMaterial.steel()

layers = m.Layer.generate_layers(thickness, [material])

N = 200
M = 4

bc = m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS
#bc = m.Model.FIXED_LEFT_RIGHT_EDGE

#corrugation_frequencies = range(2, 50, 2)
#corrugation_frequencies = [50]
corrugation_amplitudes = [0, 0.015, 0.03, 0.06, 0.1, 0.2, 0.25, 0.3]

y_per_cf = {}
x_per_cf = {}

for corrugation_amplitude in corrugation_amplitudes:

    geometry = g.General(width, curvature, corrugation_amplitude, corrugation_frequency)
    
    norm_koef = 0.5
    
#    result2D = solveLinear2D(geometry, layers, N, M, bc)
    result2D = solveLinear1D2O(geometry, layers, N, bc)
    
    
    print('========================')
    print(corrugation_amplitude)
    print(result2D.freqHz())
    print('========================')
    
    x = []
    y2D = []
    
    
    for i in range(7):
        u_max = i*norm_koef*thickness
#        result2Dnl, n = solveNonlinear2D(geometry, layers, N, M, u_max, bc)
        result2Dnl, n = solveNonlinear1D2O(geometry, layers, N, u_max, bc)
        
        d = i*norm_koef
        dy2D = result2Dnl.freqHz()/result2D.freqHz()
        
        x.append(d)
        y2D.append(dy2D)
        
        print('2D {} = {}'.format(d, dy2D))
        
        
    y_per_cf[corrugation_amplitude] = y2D
    x_per_cf[corrugation_amplitude] = x
        

folder = "./results/corrugated/"

utils.save_results(folder+"y_per_cf", y_per_cf)
utils.save_results(folder+"x_per_cf", x_per_cf)

    
    

