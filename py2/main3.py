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
		print(layer)
		layers.add(layer)
		layer_top -= layer_thickness
	return layers



width = 2
curvature = 0
thickness = 0.05

corrugation_amplitude = 0.03
corrugation_frequency = 20

N = 2
M = 1

layers_count = 1

geometry = g.CylindricalPlate(width, curvature)

layers = generate_layers(thickness, layers_count, m.Material.steel())

model = m.Model(geometry, layers, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)

mesh = mesh.Mesh.generate(width, layers, N, M, model.boundary_conditions)

print(mesh.elements)

result = s.solve(model, mesh)
#result = s.solve_nonlinearity(model, mesh)

plot.plot_strain(result, 0, width, -thickness/2, thickness/2, 0)


