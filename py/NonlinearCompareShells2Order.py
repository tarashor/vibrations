import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.shells1D.secondorder.shellsolver as s
import fem.shells1D.secondorder.result1D as r

import fem.shells1D.secondorder.nonlinearshellsolver as s_nl
import fem.shells1D.secondorder.result1Dnonlinear as r_nl


import fem.mesh as me
import plot
import numpy as np

from fem.shells1D.secondorder.matrices1D import stiffness_matrix, mass_matrix, stiffness_matrix_nl_1, stiffness_matrix_nl_2


def solveLinear(geometry, thickness, material, N, u_max, bc):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, bc)
    mesh = me.Mesh.generate2D(geometry.width, layers, N, M, model.boundary_conditions)
    
    lam, vec = s.solve(model, mesh, stiffness_matrix, mass_matrix, u_max)
    
    results_index = 10
    results = r.Result.convert_to_results(lam, vec, mesh, geometry)
    
    return results[results_index]


def solveNonlinear(geometry, thickness, material, N, M, u_max, bc):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, bc)
    mesh = me.Mesh.generate2D(geometry.width, layers, N, M, model.boundary_conditions)
    
    lam_nl, res, U1, U2, U3, um = s_nl.solve_nl(model, mesh, stiffness_matrix, mass_matrix, stiffness_matrix_nl_1, stiffness_matrix_nl_2, u_max, 10)
    
#    print('norms')
#    print(np.linalg.norm(res))
#    print(np.linalg.norm(U1))
#    print(np.linalg.norm(U2))
#    print(np.linalg.norm(U3))
#    print(np.linalg.norm(U1+U2+U3))
    
    print('maxs')
    print(get_max_u3(res, mesh))
    print(get_max_u3(U1, mesh))
    print(get_max_u3(U2, mesh))
    print(get_max_u3(U3, mesh))
    print(get_max_u3(res-U1-U2-U3, mesh))
    
    
    return r_nl.ResultNL.convert_to_result(lam_nl, res, U1, U2, U3, mesh, geometry)

def get_max_u3(v, mesh):
    u3 = v[mesh.nodes_count():2 * mesh.nodes_count()]
    Wni = np.argmax(np.absolute(u3))
    return u3[Wni]

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
thickness = 0.01

corrugation_amplitude = 0
corrugation_frequency = 0

geometry = g.General(width, curvature, corrugation_amplitude, corrugation_frequency)

N = 100
M = 4

#bc = m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS

bc = m.Model.FIXED_LEFT_RIGHT_EDGE

norm_koef = 1
u_max = norm_koef*thickness

result = solveLinear(geometry, thickness, material, N, M, u_max, bc)

resultNl = solveNonlinear(geometry, thickness, material, N, M, u_max, bc)
#resultNl2 = solveNonlinear2(geometry, thickness, material, N, M, u_max)

tN = 400
T = 6

#print(result.freqHz())
print("freqNL = {}".format(resultNl.freq))
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
