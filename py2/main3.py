import fem.geometry as g
import fem.model as m
import fem.solver as s
import fem.mesh as mesh
import plot

def generate_layers(thickness, layers_count, material):
	layer_top = thickness / 2
	layer_thickness = thickness / layers_count
	layers = set()
	for i in range(layers_count):
		layer = m.Layer(layer_top - layer_thickness, layer_top, material, i)
		layers.add(layer)
		layer_top -= layer_thickness
	return layers


#r=2
#width = r*2*3.14
#curvature = 1/r

width = 2
curvature = 0.8
thickness = 0.05

corrugation_amplitude = 0.03
corrugation_frequency = 0
#corrugation_amplitude = 0.5*thickness
#corrugation_frequency = 10

N = 50
M = 4

layers_count = 1


#geometry = g.CylindricalPlate(width, curvature)
#geometry = g.Geometry()

layers = generate_layers(thickness, layers_count, m.Material.steel())

boundary_conditions = m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS

mesh = mesh.Mesh.generate(width, layers, N, M, boundary_conditions)

corrugation_frequencies = [2, 4, 6, 8, 10, 12, 16, 20, 26, 50, 80, 100, 200, 300]

g_v = []
w_min = []

for c_f in corrugation_frequencies:

    geometry = g.CorrugatedCylindricalPlate(width, curvature, corrugation_amplitude, c_f)
    
    model = m.Model(geometry, layers, boundary_conditions)
    
    results = s.solve(model, mesh)
    print("g_v = {}, w_min = {}".format(c_f, results[0].freq))
    g_v.append(c_f)
    w_min.append(results[0].freq)
    
plt.plot(g_v, w_min, 'o-')
plt.xlabel(r"$g_v$")
plt.ylabel(r"$\omega_{min}$")
plt.title(r"Залежність $\omega_{min}$ від $g_v$" + r"($N={}, M={}$)".format(N, M))
plt.grid()
plt.show()


#result = s.solve_nonlinearity(model, mesh)


#plot.plot_strain(result, 0, width, -thickness/2, thickness/2, 0)
#results_index = 0
#plot.plot_init_geometry(geometry, results[results_index], 0, width, -thickness/2, thickness/2, 0)
#
#to_print = 20
#if (len(results) < to_print):
#    to_print = len(results)
#
#for i in range(to_print):
#    print(results[i].freq)


