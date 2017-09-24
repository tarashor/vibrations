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


def plot_displacement_norm(v1, v2, nodes, layers_count, N, M):
    x = set()
    y = set()

    list_nodes = sorted(nodes, key=lambda n: n.index)

    v = np.zeros((layers_count * M + 1, N + 1))

    for n in list_nodes:
        x.add(n.x)
        y.add(n.y)
        i = n.index // (N + 1)
        j = n.index % (N + 1)
        norm = np.sqrt(v1[n.index] * v1[n.index] + v2[n.index] * v2[n.index])
        v[i, j] = norm

    x = sorted(x)
    y = sorted(y)

    (X, Y) = np.meshgrid(x, y)
    surf = plt.contourf(X, Y, v, cmap=cm.rainbow)
    plt.colorbar(surf)
    plt.show()


def plot_init_geometry(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, N, M):
    l, v1, v2, nodes = get_lowest_freq(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, N, M)

    X_bottom = []
    Y_bottom = []
    X_top = []
    Y_top = []

    lnodes = sorted(nodes, key=lambda n: n.index)

    for i in range(N + 1):
        ind = (layers_count * M) * (N + 1) + i
        # ind = i
        a1 = lnodes[ind].x
        a2 = lnodes[ind].y
        ar = (np.pi + curvature * width) / 2 - a1 * curvature
        x = (1 / curvature + a2 + corrugation_amplitude * np.cos(corrugation_frequency * ar)) * np.cos(ar)
        y = (1 / curvature + a2 + corrugation_amplitude * np.cos(corrugation_frequency * ar)) * np.sin(ar)
        X_bottom.append(x)
        Y_bottom.append(y)
        a1 += v2[lnodes[ind].index]
        a2 += v1[lnodes[ind].index]
        ar = (np.pi + curvature * width) / 2 - a1 * curvature
        x = (1 / curvature + a2 + corrugation_amplitude * np.cos(corrugation_frequency * ar)) * np.cos(ar)
        y = (1 / curvature + a2 + corrugation_amplitude * np.cos(corrugation_frequency * ar)) * np.sin(ar)
        X_top.append(x)
        Y_top.append(y)

    plt.plot(X_bottom, Y_bottom, label="UNdeformed")
    plt.plot(X_top, Y_top, label="deformed")
    plt.legend(loc='best')
    plt.show()


def plot_sample(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, N, M):
    l, v1, v2, nodes = get_lowest_freq(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, N, M)
    print(l)
    plot_displacement_norm(v1, v2, nodes, layers_count, N, M)


def calculate_data_freq_from_NxM(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, N_max, M_max):
    data = []

    for n in range(10, N_max + 1, 10):
        for m in range(2, M_max + 1, 2):
            l, v1, v2, nodes = get_lowest_freq(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, n, m)
            data.append([n, m, l])

    return data


def calculate_data_freq_from_layers_count(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, N, M, layers_count_max):
    data = []

    for lc in range(layers_count_max + 1):
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

    ax.set_yticks(np.arange(2, 21, 2))
    ax.set_xlabel(r'$N$', fontsize=14)
    ax.set_ylabel(r'$M$', fontsize=14)
    ax.set_zlabel(r'$\omega_{min}  $', fontsize=14)
    ax.zaxis.set_rotate_label(False)

    ax.title.set_text(r'Залежність $\omega_{min}$ від к-сті елементів по товщині $M$ і по довжині $N$ (к-сть шарів = 1)')

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()


def plot_freq_from_layers_count(data):
    pass


width = 2
curvature = 0.8
thickness = 0.05

corrugation_amplitude = 0
corrugation_frequency = 0


freq_from_NM_file = "freq_from_NxM"
freq_from_layers_file = "freq_from_layers_count"

freq_from_NM_file_done = "freq_from_NxM_done"
freq_from_layers_file_done = "freq_from_layers_count_done"

layers_count_default = 1
N_default = 20
M_default = 2

# 1
#data = calculate_data_freq_from_NxM(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count_default, 100, 20)
#utils.save_in_file(freq_from_NM_file, data)

# 2
#data = utils.read_from_file(freq_from_NM_file_done)
# plot_freq_from_NxM(data)

# 3
#data = calculate_data_freq_from_layers_count(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, N_default, M_default, 5)
#utils.save_in_file(freq_from_layers_file, data)

# 4
#data = utils.read_from_file(freq_from_layers_file_done)
# plot_freq_from_layers_count(data)

# 5
#plot_sample(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count_default, N_default, M_default)

# 6
plot_init_geometry(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count_default, N_default, M_default)
