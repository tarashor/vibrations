import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

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