import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.general2D.solver_nonlinear as s_nl
import fem.general2D.result2Dnonlinear as r
import fem.mesh as me
import plot

from fem.general2D.matrices2D import stiffness_matrix, mass_matrix, stiffness_matrix_nl_1, stiffness_matrix_nl_2


def solve(geometry, thickness, material, N, M, u_max):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
#    model = m.Model(geometry, layers, m.Model.FIXED_LEFT_RIGHT_EDGE)
    mesh = me.Mesh.generate2D(geometry.width, layers, N, M, model.boundary_conditions)
    
    lam_nl, res, U1, U2, U3 = s_nl.solve_nl(model, mesh, stiffness_matrix, mass_matrix, stiffness_matrix_nl_1, stiffness_matrix_nl_2, u_max)
    
    return mesh, lam_nl, res, U1, U2, U3

    
#E = 40000
#v = 0.3
#rho = 2000
#
#material = mat.IsotropicMaterial(E,v,rho)
    
material = mat.IsotropicMaterial.steel()

width = 1
curvature = 0
thickness = 0.1

corrugation_amplitude = 0
corrugation_frequency = 0

geometry = g.General(width, curvature, corrugation_amplitude, corrugation_frequency)

N = 100
M = 4

norm_koef = 0.1
u_max = norm_koef*thickness

mesh, lam_nl, res, U1, U2, U3 = solve(geometry, thickness, material, N, M, u_max)
result = r.ResultNL.convert_to_result(lam_nl, res, U1, U2, U3, mesh, geometry)

tN = 1000
T = 6

#plot.plot_animate_in_cartesian(result, 0, width, -thickness / 2, thickness / 2, T, tN, geometry.to_cartesian_coordinates)
#
##for t in range (0,tN,100):
##    time = T/tN*t
##    plot.plot_init_and_deformed_geometry_in_cartesian(result, 0, width, -thickness / 2, thickness / 2, time, geometry.to_cartesian_coordinates)

print(result.freqHz())


x = []
y = []

deltat = T / tN

for t in range(tN):
    time = deltat*t
    u = result.get_displacement_and_deriv(width/2, 0, 0, time)
    u1 = u[0]
    u2 = u[4]
    u3 = u[8]
    x.append(time)
    y.append(u3)
    

plot.plot_vibrations(x,y)
