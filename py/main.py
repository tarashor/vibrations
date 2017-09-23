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


def get_lowest_freq(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, N, M):

    layer_top = thickness / 2
    layer_thickness = thickness / layers_count
    layers = set()
    for i in range(layers_count):
        layer = fem.model.Layer(layer_top - layer_thickness, layer_top, fem.model.Material.steel(), i)
        layers.add(layer)
        layer_top -= layer_thickness

    result = solve(width, curvature, corrugation_amplitude, corrugation_frequency, layers, N, M)

    return result.get_result(0)

def plot_sample(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, N, M):
    l, v1, v2, nodes = get_lowest_freq(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, N, M)
    print(l)

    x = []
    y = []

    list_nodes = sorted(nodes, key=lambda n: n.index)

    for n in list_nodes:
        x.append(n.x + v1[n.index])
        y.append(n.y + v2[n.index])

    plt.plot(x, y, 'ro')
    plt.show()

def calculate_data_freq_from_NxM(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, N_max, M_max):
    data = []

    for n in range(10, N_max+1, 10):
        for m in range(2, M_max+1, 2):
            l, v1, v2, nodes = get_lowest_freq(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, n, m)
            data.append([n,m,l])
            
    return data

def calculate_data_freq_from_layers_count(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, N, M, layers_count_max):
    data = []

    for lc in range(layers_count_max+1):
        l, v1, v2, nodes = get_lowest_freq(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, lc, N, M)
        data.append([lc, l])
            
    return data
    


def plot_freq_from_NxM(data):
    # Make data.
    n = set()
    m = set()
    for line in data:
        n.add(int(line[0]))
        m.add(int(line[1]))

    n = sorted(n)
    m = sorted(m)

    z = np.zeros((len(n), len(m)))
    for i in range(len(n)):
        for j in range(len(m)):
            freq = 0
            for item in data:
                if (int(item[0]) == n[i] and int(item[1]) == m[j]):
                    freq = float(item[2])
                    break
            z[i, j] = freq

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    n, m = np.meshgrid(n, m)
    surf = ax.plot_surface(n, m, z, cmap=cm.rainbow)

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()
    
    
def plot_freq_from_layers_count(data):
    pass


width = 1
curvature = 0.25
thickness = 0.08

corrugation_amplitude = 0.03
corrugation_frequency = 20


freq_from_NM_file = "freq_from_NxM"
freq_from_layers_file = "freq_from_layers_count"

freq_from_NM_file_done = "freq_from_NxM_done"
freq_from_layers_file_done = "freq_from_layers_count_done"

layers_count_default = 1
N_default = 70
M_default = 10

# 1
#data = calculate_data_freq_from_NxM(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count_default, 100, 20)
#utils.save_in_file(freq_from_NM_file, data)

# 2
data = utils.read_from_file(freq_from_NM_file_done)
plot_freq_from_NxM(data)

# 3
#data = calculate_data_freq_from_layers_count(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, N_default, M_default, 5)
#utils.save_in_file(freq_from_layers_file, data)

# 4
#data = utils.read_from_file(freq_from_layers_file_done)
#plot_freq_from_layers_count(data)

# 5
#plot_sample(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count_default, N_default, M_default)
