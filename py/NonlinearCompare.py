import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.general2D.solverlinear as s
import fem.general2D.result2D as r

import fem.general2D.solver_nonlinear as s_nl
import fem.general2D.result2Dnonlinear as r_nl

import fem.general2D.solver_nonlinear2 as s_nl2

import fem.mesh as me
import plot

from fem.general2D.matrices2D import stiffness_matrix, mass_matrix, stiffness_matrix_nl_1, stiffness_matrix_nl_2


def solveLinear(geometry, thickness, material, N, M, u_max, bc):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, bc)
    mesh = me.Mesh.generate2D(geometry.width, layers, N, M, model.boundary_conditions)
    
    lam, vec = s.solve(model, mesh, stiffness_matrix, mass_matrix, u_max)
    
    results_index = 0
    results = r.Result.convert_to_results(lam, vec, mesh, geometry)
    
    return results[results_index]


def solveNonlinear(geometry, thickness, material, N, M, u_max, bc):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, bc)
    mesh = me.Mesh.generate2D(geometry.width, layers, N, M, model.boundary_conditions)
    
    lam_nl, res, U1, U2, U3, um = s_nl.solve_nl(model, mesh, stiffness_matrix, mass_matrix, stiffness_matrix_nl_1, stiffness_matrix_nl_2, u_max)
    
    return r_nl.ResultNL.convert_to_result(lam_nl, res, U1, U2, U3, mesh, geometry)

#def solveNonlinear2(geometry, thickness, material, N, M, u_max):
#    layers = m.Layer.generate_layers(thickness, [material])
#    model = m.Model(geometry, layers, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
##    model = m.Model(geometry, layers, m.Model.FIXED_LEFT_RIGHT_EDGE)
#    mesh = me.Mesh.generate2D(geometry.width, layers, N, M, model.boundary_conditions)
#    
#    lam_nl, res = s_nl2.solve_nl(model, mesh, stiffness_matrix, mass_matrix, stiffness_matrix_nl_1, stiffness_matrix_nl_2, u_max)
#    
#    return r.Result.convert_to_result(lam_nl, res, mesh, geometry)


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

#bc = m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS

bc = m.Model.FIXED_LEFT_RIGHT_EDGE

norm_koef = 2
u_max = norm_koef*thickness

result = solveLinear(geometry, thickness, material, N, M, u_max, bc)

resultNl = solveNonlinear(geometry, thickness, material, N, M, u_max, bc)
#resultNl2 = solveNonlinear2(geometry, thickness, material, N, M, u_max)

tN = 400
T = 6

print(result.freqHz())
print(resultNl.freqHz())
#print(resultNl2.freqHz())


x = []
y = []
ynl = []
ynl2 = []

deltat = T / tN

for t in range(tN):
    time = deltat*t
    u = result.get_displacement_and_deriv(width/2, 0, 0, time)
    unl = resultNl.get_displacement_and_deriv(width/2, 0, 0, time)
#    unl2 = resultNl2.get_displacement_and_deriv(width/2, 0, 0, time)
    x.append(time)
    y.append(u[8])
    ynl.append(unl[8])
#    ynl2.append(unl2[8])
    
    

plot.plot_vibrations(x,y,ynl)
