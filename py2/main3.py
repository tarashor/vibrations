import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.solver as s
import fem.mesh as me
import plot

from fem.matrices import stiffness_matrix, mass_matrix


def generate_layers(thickness, layers_count, material):
    layer_top = thickness / 2
    layer_thickness = thickness / layers_count
    layers = set()
    for i in range(layers_count):
        layer = m.Layer(layer_top - layer_thickness, layer_top, material, i)
        layers.add(layer)
        layer_top -= layer_thickness
    return layers


def solve(width, curvature, thickness, corrugation_amplitude, corrugation_frequency):
    layers_count = 1
    layers = generate_layers(thickness, layers_count, mat.IsotropicMaterial.steel())
    mesh = me.Mesh.generate(width, layers, N, M, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
#    geometry = g.CorrugatedCylindricalPlate(width, curvature, corrugation_amplitude, corrugation_frequency)
    geometry = g.CylindricalPlate(width, curvature)
#    geometry = g.Geometry()
    model = m.Model(geometry, layers, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
    return s.solve(model, mesh, stiffness_matrix, mass_matrix)


# r=2
# width = r*2*3.14
# curvature = 1/r

width = 2
curvature = 0.8
thickness = 0.05

corrugation_amplitude = 0.03
corrugation_frequency = 20
# corrugation_amplitude = 0.5*thickness
# corrugation_frequency = 10

N = 50
M = 4


results = solve(width, curvature, thickness, corrugation_amplitude, corrugation_frequency)
results_index = 0
#plot.plot_mesh(results[results_index].mesh, width, thickness)

#print(results[results_index].u1)

#plot.plot_deformed_mesh(results[results_index], width, thickness)

plot.plot_init_and_deformed_geometry(results[results_index], 0, width, -thickness / 2, thickness / 2, 0)
#plot.plot_init_geometry(results[results_index].geometry, 0, width, -thickness / 2, thickness / 2, 0)
# plot.plot_strain(results[results_index], 0, width, -thickness / 2, thickness / 2, 0)


to_print = 20
if (len(results) < to_print):
    to_print = len(results)

for i in range(to_print):
    print(results[i].freq)
