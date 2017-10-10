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

        # v[i, j] = v2[n.index]

    x = sorted(x)
    y = sorted(y)

    (X, Y) = np.meshgrid(x, y)
    surf = plt.contourf(X, Y, v, cmap=cm.rainbow)
    plt.colorbar(surf)
    plt.show()


def plot_init_geometry(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, N, M):
    l, v1, v2, nodes = get_lowest_freq(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, N, M)

    print("Min freq = {}".format(l))
    print("v1 = {}".format(v1))
    print("v2 = {}".format(v2))

    lnodes = sorted(nodes, key=lambda n: n.index)
    # d1 = []
    # d2 = []
    # for n in lnodes:
    #     d1.append(n.x + v2[n.index])
    #     d2.append(n.y + v1[n.index])

    # plt.plot(d1, d2, 'ro')
    # # plt.axis([-2, 2, -2, 2])
    # plt.show()

    X_init = []
    Y_init = []
    X_deformed = []
    Y_deformed = []

    for i in range(N + 1):
        ind = (layers_count * M) * (N + 1) + i
        # ind = i

        x = lnodes[ind].x
        y = lnodes[ind].y

        if (curvature > 0):
            ar = (np.pi + curvature * width) / 2 - x * curvature
            x = (1 / curvature + y + corrugation_amplitude * np.cos(corrugation_frequency * ar)) * np.cos(ar)
            y = (1 / curvature + y + corrugation_amplitude * np.cos(corrugation_frequency * ar)) * np.sin(ar)

        X_init.append(x)
        Y_init.append(y)

        x = lnodes[ind].x + v1[lnodes[ind].index]
        y = lnodes[ind].y + v2[lnodes[ind].index]
        if (curvature > 0):
            ar = (np.pi + curvature * width) / 2 - x * curvature
            x = (1 / curvature + y + corrugation_amplitude * np.cos(corrugation_frequency * ar)) * np.cos(ar)
            y = (1 / curvature + y + corrugation_amplitude * np.cos(corrugation_frequency * ar)) * np.sin(ar)

        X_deformed.append(x)
        Y_deformed.append(y)

    for i in range(N + 1):
        ind = N - i

        x = lnodes[ind].x
        y = lnodes[ind].y

        if (curvature > 0):
            ar = (np.pi + curvature * width) / 2 - x * curvature
            x = (1 / curvature + y + corrugation_amplitude * np.cos(corrugation_frequency * ar)) * np.cos(ar)
            y = (1 / curvature + y + corrugation_amplitude * np.cos(corrugation_frequency * ar)) * np.sin(ar)

        X_init.append(x)
        Y_init.append(y)

        x = lnodes[ind].x + v1[lnodes[ind].index]
        y = lnodes[ind].y + v2[lnodes[ind].index]
        if (curvature > 0):
            ar = (np.pi + curvature * width) / 2 - x * curvature
            x = (1 / curvature + y + corrugation_amplitude * np.cos(corrugation_frequency * ar)) * np.cos(ar)
            y = (1 / curvature + y + corrugation_amplitude * np.cos(corrugation_frequency * ar)) * np.sin(ar)

        X_deformed.append(x)
        Y_deformed.append(y)

    X_init.append(X_init[0])
    Y_init.append(Y_init[0])
    X_deformed.append(X_deformed[0])
    Y_deformed.append(Y_deformed[0])

    plt.plot(X_init, Y_init, label="початкова конфігурація")
    plt.plot(X_deformed, Y_deformed, label="поточна конфігурація")
    plt.title("Деформації при {}-ій власній частоті".format(3))
    # plt.title(r"Форма панелі з такими параметрами $l={}, h={}, K={}, g_A={}, g_v={}$".format(width, thickness, curvature, corrugation_amplitude, corrugation_frequency))
    # plt.axis([-1, 1, 3, 5])
    plt.legend(loc='best')
    plt.grid()
    plt.show()


def plot_sample(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, N, M):
    l, v1, v2, nodes = get_lowest_freq(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, N, M)
    print(l)
    plot_displacement_norm(v1, v2, nodes, layers_count, N, M)


def plot_freq_from_corrugation_frequency(width, thickness, curvature, corrugation_amplitude, corrugation_frequencies, layers_count, N, M):
    cfs = []
    ls = []

    for cf in corrugation_frequencies:
        l, v1, v2, nodes = get_lowest_freq(width, thickness, curvature, corrugation_amplitude, cf, layers_count, N, M)
#        print("{},{}".format(cf, l))
#        data.append([cf, l])

        cfs.append(cf)
        ls.append(l)

    plt.plot(cfs, ls, 'o-')
    plt.xlabel(r"$g_v$")
    plt.ylabel(r"$\omega_{min}$")
    plt.title(r"Залежність $\omega_{min}$ від $g_v$" + r"($N={}, M={}$)".format(N, M))
    plt.grid()
    plt.show()


def calculate_data_freq_from_NxM(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, N_max, M_max):
    data = []

    for n in range(40, N_max + 1, 40):
        for m in range(4, M_max + 1, 4):
            l, v1, v2, nodes = get_lowest_freq(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count, n, m)
            print("{},{},{}".format(n, m, l))
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
    print(data)
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
    m, n = np.meshgrid(m, n)

    print(n)
    print(m)
    print(z)

    surf = ax.plot_surface(n, m, z, cmap=cm.rainbow)

    ax.set_yticks(np.arange(4, 17, 2))
    ax.set_xlabel(r'$N$', fontsize=14)
    ax.set_ylabel(r'$M$', fontsize=14)
    ax.set_zlabel(r'$\omega_{min}  $', fontsize=14)
    ax.zaxis.set_rotate_label(False)

    ax.title.set_text(r'Залежність $\omega_{min}$ від к-сті елементів по товщині $M$ і по довжині $N$ (к-сть шарів = 1)')

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()


def plot_freq_from_layers_count(data):
    x = []
    y = []
    for line in data:
        x.append(int(line[0]))
        y.append(float(line[1]))

    plt.plot(x, y)
    plt.xlabel("Кількість шарів")
    plt.ylabel(r"$\omega_{min}$")
    plt.title(r"Залежність $\omega_{min}$ від к-сті шарів сталі ($M=70, N=10, h=0.05$)")
    plt.grid()
    plt.show()


width = 2
curvature = 0
thickness = 0.05

corrugation_amplitude = 0.03
corrugation_frequency = 20


freq_from_NM_file = "freq_from_NxM"
freq_from_layers_file = "freq_from_layers_count"

freq_from_NM_file_done = "freq_from_NxM_done"
freq_from_layers_file_done = "freq_from_layers_count_done"

layers_count_default = 1
N_default = 20
M_default = 4

# 1
# data = calculate_data_freq_from_NxM(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count_default, 240, 16)
# utils.save_in_file(freq_from_NM_file, data)

# 2
# data = utils.read_from_file(freq_from_NM_file_done)
# plot_freq_from_NxM(data)

# 3
# data = calculate_data_freq_from_layers_count(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, N_default, M_default, 5)
# utils.save_in_file(freq_from_layers_file, data)

# 4
# data = utils.read_from_file(freq_from_layers_file_done)
# plot_freq_from_layers_count(data)

# 5
# plot_sample(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count_default, N_default, M_default)

# 6
plot_init_geometry(width, thickness, curvature, corrugation_amplitude, corrugation_frequency, layers_count_default, N_default, M_default)

# 7
# plot_freq_from_corrugation_frequency(width, thickness, curvature, corrugation_amplitude, [2,4,6,8,10,12,16,20,26,50,80,100], layers_count_default, N_default, M_default)
