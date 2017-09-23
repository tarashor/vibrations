import fem.model
import fem.mesh
import fem.solver
import csv
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


def get_lowest_freq(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, N, M):
    geometry = fem.model.Geometry(width, curvature, corrugation_amplitude, corrugation_frequency)

    layer_top = thickness / 2
    layer_thickness = thickness / layers_count
    layers = set()
    for i in range(layers_count):
        layer = fem.model.Layer(layer_top - layer_thickness, layer_top, fem.model.Material.steel(), i)
        layers.add(layer)
        layer_top -= layer_thickness

    model = fem.model.Model(geometry, layers, fem.model.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)

    mesh = fem.mesh.Mesh.generate(model.geometry.width, layers, N, M, model.boundary_conditions)

    result = fem.solver.solve(model, mesh)

    l, v1, v2 = result.get_result(0)
    return l


def save_in_file(file_name, data):
    with open(file_name, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for line in data:
            writer.writerow(line)


def read_from_file(file_name):
    data = []
    with open(file_name, 'r') as csvfile:
        reader = csv.reader(csvfile)
        for line in reader:
            data.append(line)

    return data


def plot_converg(data):
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


width = 1
curvature = 0.25
thickness = 0.08

corrugation_amplitude = 0.03
corrugation_frequency = 20

layers_count = 1

N = 70
M = 10

freq_from_NM_file = "data.csv"

# data = []
#
# for n in range(10, N+1, 10):
#    for m in range(2, M+1, 2):
#        l = get_lowest_freq(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, n, m)
#        data.append([n,m,l])
#
# save_in_file(freq_from_NM_file, data)

# data = read_from_file(freq_from_NM_file)

# plot_converg(data)

freq_from_layers_file = "freq_from_layers_file2.csv"

data = []

for lc in range(5, 11):
    l = get_lowest_freq(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, lc, N, M)
    print("{} - {}".format(lc, l))

    data.append([lc, l])

save_in_file(freq_from_layers_file, data)


# print(l)
#
#x = []
#y = []
#
#list_nodes = sorted(mesh.nodes, key=lambda n: n.index)
#print("Nodes = {}".format(len(mesh.nodes)))
#
#fn = mesh.fixed_nodes
# for n in fn:
#    print(n.x + v1[n.index])
#    print(n.y + v2[n.index])
#
# for n in list_nodes:
#    x.append(n.x + v1[n.index])
#    y.append(n.y + v2[n.index])
#
#plt.plot(x, y, 'ro')
## plt.axis([-2, 2, -2, 2])
# plt.show()
