import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.mesh as me
import plot

import fem.shells1D.secondorder.shellsolver as s
import fem.shells1D.secondorder.result1D as r
from fem.shells1D.secondorder.matrices1D import stiffness_matrix, mass_matrix


def solve(geometry, thickness, material, N):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
    model.fixed_x3 = -thickness/2
    mesh = me.Mesh.generate1D(geometry.width, layers, N, model.boundary_conditions)
    
    lam, vec = s.solve(model, mesh, stiffness_matrix, mass_matrix)
    
    return lam, vec, mesh, geometry


    
#E = 40000000000
#v = 0.3
#rho = 2000
#
#material = mat.IsotropicMaterial(E,v,rho)
    
material = mat.IsotropicMaterial.steel()

width = 1
curvature = 2
thickness = 0.01

corrugation_amplitude = 0
corrugation_frequency = 0

geometry = g.General(width, curvature, corrugation_amplitude, corrugation_frequency)

N = 200


lam, vec, mesh, geometry = solve(geometry, thickness, material, N)
results = r.Result.convert_to_results(lam, vec, mesh, geometry, thickness)

results_index = 0

#plot.plot_mesh(results[results_index].mesh, width, thickness)


plot.plot_init_and_deformed_geometry_in_cartesian(results[results_index], 0, width, -thickness / 2, thickness / 2, 0, geometry.to_cartesian_coordinates)

to_print = 20
if (len(results) < to_print):
    to_print = len(results)

for i in range(to_print):
    print(results[i].freqHz())
