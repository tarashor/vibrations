import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.solver as s
import fem.dynamic_solver as ds
import fem.mesh as me
import numpy as np
import matplotlib.pyplot as plt
import plot


from fem.dmatrices import tangent_stiffness_matrix, mass_matrix, force_vector


def generate_layers(thickness, layers_count, material):
    layer_top = thickness / 2
    layer_thickness = thickness / layers_count
    layers = set()
    for i in range(layers_count):
        layer = m.Layer(layer_top - layer_thickness, layer_top, material, i)
        layers.add(layer)
        layer_top -= layer_thickness
    return layers


def solve_dynamic(geometry, thickness, T, time_intervals, u_1_0, u_3_0, v_1_0, v_3_0):
    layers_count = 1
    layers = generate_layers(thickness, layers_count, mat.IsotropicMaterial.steel())
    model = m.Model(geometry, layers, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
    mesh = me.Mesh.generate(width, layers, N, M, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
    nodes_count = mesh.nodes_count()
    u0 = np.zeros((2*nodes_count))
    v0 = np.zeros((2*nodes_count))
    for node in mesh.nodes:
        u0[node.index]=u_1_0(node.x1, node.x2)
        u0[node.index+nodes_count]=u_3_0(node.x1, node.x2)
        v0[node.index]=v_1_0(node.x1, node.x2)
        v0[node.index+nodes_count]=v_3_0(node.x1, node.x2)
        
    return ds.solve(model, mesh, tangent_stiffness_matrix, mass_matrix, force_vector, T, time_intervals, u0, v0)

def u_1_0(x1, x3):
    return np.sin(np.pi*x1/width)

def u_3_0(x1, x3):
    return np.sin(2*np.pi*x1/width)

def v_1_0(x1, x3):
    return 0

def v_3_0(x1, x3):
    return 0


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
time_intervals = 50

dresults = solve_dynamic(geometry, thickness, T, time_intervals, u_1_0, u_3_0, v_1_0, v_3_0)

nodes_count = (N+1)*(M+1)
u1=[]
u3=[]

for t in range(time_intervals):
    u=dresults[t]
    r1=u[nodes_count//2]
    r3=u[nodes_count//2+nodes_count]
    print(r1)
    print(r3)
    u1.append(r1)
    u3.append(r3)

time_points = np.linspace(0, T, time_intervals)

plt.plot(range(time_intervals), u1, "r", label="u1")
plt.plot(range(time_intervals), u3, "b", label="u3")
plt.title("Переміщення точки")
# plt.title(r"Форма панелі $L={}, h={}, K={}, g_A={}, g_v={}$".format(x1_end - x1_start, x2_end - x2_start, result.geometry.curvature, result.geometry.corrugation_amplitude, result.geometry.corrugation_frequency))
plt.axes().set_aspect('equal', 'datalim')
plt.legend(loc='best')
plt.xlabel(r"$час$, sec", fontsize=12)
plt.ylabel(r"$переміщення$, м", fontsize=12)
plt.grid()
plt.show()
#x1 = width/2
#x2 = 0
#x3 = 0
#time_points = np.linspace(0, T, time_intervals)
#
#plot.plot_point_in_time(dresult, x1, x2, x3, time_points)
    

    

    