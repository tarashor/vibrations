import fem.model
import fem.mesh
import fem.solver
import utils
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


def solve(width, curvature, corrugation_amplitude, corrugation_frequency, layers, N, M):
    geometry = fem.model.Geometry(width, curvature, corrugation_amplitude, corrugation_frequency)

    model = fem.model.Model(geometry, layers, fem.model.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)

    mesh = fem.mesh.Mesh.generate(model.geometry.width, layers, N, M, model.boundary_conditions)

    return fem.solver.solve(model, mesh)
#    return fem.solver.solve_nonlinearity(model, mesh)


def get_result_for_same_layers(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, N, M):

    layer_top = thickness / 2
    layer_thickness = thickness / layers_count
    layers = set()
    for i in range(layers_count):
        layer = fem.model.Layer(layer_top - layer_thickness, layer_top, fem.model.Material.steel(), i)
        layers.add(layer)
        layer_top -= layer_thickness

    return solve(width, curvature, corrugation_amplitude, corrugation_frequency, layers, N, M)


    
def plot_grad_norm(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, N, M):
    result = get_result_for_same_layers(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, N, M)

    freq_index = 0
    x = set()
    y = set()

    list_nodes = sorted(result.get_nodes(), key=lambda n: n.index)

    v = np.zeros((layers_count * M + 1, N + 1))

    for n in list_nodes:
        x.add(n.x)
        y.add(n.y)
        i = n.index // (N + 1)
        j = n.index % (N + 1)
        v[i, j] = result.get_gradu(freq_index, n.x, n.y)[0]

        # v[i, j] = v2[n.index]

    x = sorted(x)
    y = sorted(y)

    (X, Y) = np.meshgrid(x, y)
    surf = plt.contourf(X, Y, v, cmap=cm.rainbow)
    plt.colorbar(surf)
    plt.show()



    
 # -*- coding: utf-8 -*-

width = 2
curvature = 0
thickness = 0.05

corrugation_amplitude = 0.03
corrugation_frequency = 20

layers_count_default = 1
N_default = 40
M_default = 4

plot_init_geometry(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count_default, N_default, M_default)