import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from fem import geometry as g

plot_x1_elements = 400
plot_x2_elements = 20


def plot_strain(result, x1_start, x1_end, x2_start, x2_end, time):

    dx1 = (x1_end - x1_start) / plot_x1_elements
    dx2 = (x2_end - x2_start) / plot_x2_elements
    
    x1 = set()
    x2 = set()
    
    v = np.zeros((plot_x1_elements+1, plot_x2_elements+1))
    
    for i in range(plot_x1_elements + 1):
        for j in range(plot_x2_elements + 1):
            x1c = x1_start + i * dx1
            x2c = x2_start + j * dx2
            v[i, j] = result.get_strain(x1c, x2c, 0, time)[3]
            x1.add(x1c)
            x2.add(x2c)

    x1 = sorted(x1)
    x2 = sorted(x2)

    (X1, X2) = np.meshgrid(x1, x2)
    surf = plt.contourf(X1, X2, v.T, cmap=cm.rainbow)
    plt.colorbar(surf)
    plt.show()
    
def plot_init_geometry(geometry, result, x1_start, x1_end, x2_start, x2_end, time):

    dx1 = (x1_end - x1_start) / plot_x1_elements
    
    X_init = []
    Y_init = []
    X_deformed = []
    Y_deformed = []

    x3 = 0
    x2 = x2_end
    for i in range(plot_x1_elements + 1):
        x1 = x1_start + i * dx1
        u = result.get_displacement_and_deriv(x1, x2, x3, time)
        u1 = u[0]
        u2 = u[4]
        
        x, y, z = geometry.to_cartesian_coordinates(x1, x2, x3)

        X_init.append(x)
        Y_init.append(y)

        x1_deformed = x1 + u1
        x2_deformed = x2 + u2
        
        x, y, z = geometry.to_cartesian_coordinates(x1_deformed, x2_deformed, x3)

        X_deformed.append(x)
        Y_deformed.append(y)

    x2 = x2_start
    for i in range(plot_x1_elements + 1):
        x1 = x1_end - i * dx1
        u = result.get_displacement_and_deriv(x1, x2, x3, time)
        u1 = u[0]
        u2 = u[4]
        
        x, y, z = geometry.to_cartesian_coordinates(x1, x2, x3)

        X_init.append(x)
        Y_init.append(y)

        x1_deformed = x1 + u1
        x2_deformed = x2 + u2
        
        x, y, z = geometry.to_cartesian_coordinates(x1_deformed, x2_deformed, x3)

        X_deformed.append(x)
        Y_deformed.append(y)

    X_init.append(X_init[0])
    Y_init.append(Y_init[0])
    X_deformed.append(X_deformed[0])
    Y_deformed.append(Y_deformed[0])

    plt.plot(X_init, Y_init, label="початкова конфігурація")
    plt.plot(X_deformed, Y_deformed, label="поточна конфігурація")
    plt.title("Переміщення")
    # plt.title(r"Форма панелі з такими параметрами $l={}, h={}, K={}, g_A={}, g_v={}$".format(width, thickness, curvature, corrugation_amplitude, corrugation_frequency))
    # plt.axis([-1, 1, 3, 5])
    plt.legend(loc='best')
    plt.grid()
    plt.show()