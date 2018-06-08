import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.general2D.solverlinear as s
import fem.general2D.result2D as r
import fem.mesh as me
import plot

from fem.general2D.matrices2D import stiffness_matrix, mass_matrix, stiffness_matrix_nl


def solve(geometry, thickness, material, N, M):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
    mesh = me.Mesh.generate2D(geometry.width, layers, N, M, model.boundary_conditions)
    
    
    lam, vec = s.solve(model, mesh, stiffness_matrix, mass_matrix)
    
    return lam, vec, mesh, geometry


    
material = mat.IsotropicMaterial.steel()

width = 2
curvature = 0.000000001
thickness = 0.05

corrugation_amplitude = 0
corrugation_frequency = 20

geometry = g.General(width, curvature, corrugation_amplitude, corrugation_frequency)

N = 100
M = 4


lam, vec, mesh, geometry = solve(geometry, thickness, material, N, M)
results = r.Result.convert_to_results(lam, vec, mesh, geometry)

results_index = 0


plot.plot_init_and_deformed_geometry_in_cartesian(results[results_index], 0, width, -thickness / 2, thickness / 2, 0, geometry.to_cartesian_coordinates)

to_print = 20
if (len(results) < to_print):
    to_print = len(results)

for i in range(to_print):
    print(results[i].freqHz())
