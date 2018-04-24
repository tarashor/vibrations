import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.solver as s
import fem.dynamic_solver as ds
import fem.mesh as me
import numpy as np
import plot


from fem.matrices import stiffness_matrix, mass_matrix, stiffness_matrix_nl


def generate_layers(thickness, layers_count, material):
    layer_top = thickness / 2
    layer_thickness = thickness / layers_count
    layers = set()
    for i in range(layers_count):
        layer = m.Layer(layer_top - layer_thickness, layer_top, material, i)
        layers.add(layer)
        layer_top -= layer_thickness
    return layers


def solve_dynamic(geometry, thickness, T, time_intervals):
    layers_count = 1
    layers = generate_layers(thickness, layers_count, mat.IsotropicMaterial.steel())
    model = m.Model(geometry, layers, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
    mesh = me.Mesh.generate(width, layers, N, M, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
    return ds.solve(model, mesh, stiffness_matrix, mass_matrix, stiffness_matrix_nl, T, time_intervals)


# r=2
# width = r*2*3.14
# curvature = 1/r

width = 2
curvature = 1
thickness = 0.05

corrugation_amplitude = 0.03
corrugation_frequency = 20

#geometry = g.CorrugatedCylindricalPlate(width, curvature, corrugation_amplitude, corrugation_frequency)
geometry = g.CylindricalPlate(width, curvature)
#geometry = g.Plate()

N = 50
M = 4

T = 1
time_intervals = 200

dresult = solve_dynamic(geometry, thickness, T, time_intervals)

x1 = width/2
x2 = 0
x3 = 0
time_points = np.linspace(0, T, time_intervals)

plot.plot_point_in_time(dresult, x1, x2, x3, time_points)
    

    

    