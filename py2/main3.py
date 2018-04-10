import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.solver as s
import fem.mesh as me
import plot

import pickle

from fem.matrices import stiffness_matrix, mass_matrix


def generate_layers(thickness, layers_count, material):
    layer_top = thickness / 2
    layer_thickness = thickness / layers_count
    layers = set()
    for i in range(layers_count):
        layer = m.Layer(layer_top - layer_thickness, layer_top, material, i)
        layers.add(layer)
        layer_top -= layer_thickness
    return layers


def solve(geometry, thickness):
    layers_count = 1
    layers = generate_layers(thickness, layers_count, mat.IsotropicMaterial.steel())
    model = m.Model(geometry, layers, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
    mesh = me.Mesh.generate(width, layers, N, M, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
    lam, vec = s.solve(model, mesh, stiffness_matrix, mass_matrix)
    
    return lam, vec, mesh, geometry


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

toCalculate = False

def save_mesh(filename, mesh):
    with open(filename + '.mesh', 'wb') as f:
        pickle.dump(mesh, f)
        
def save_geometry(filename, geometry):
    with open(filename + '.geom', 'wb') as f:
        pickle.dump(geometry, f)

def load_geometry(filename):
    with open(filename + '.geom', 'rb') as f:
        return pickle.load(f)
    
def load_mesh(filename):
    with open(filename + '.mesh', 'rb') as f:
        return pickle.load(f)
    
def save_results(filename, results):
    with open(filename + '.res', 'wb') as f:
        pickle.dump(results, f)

def load_results(filename):
    with open(filename + '.res', 'rb') as f:
        return pickle.load(f)
    
filename = str(geometry) + "_{}x{}".format(N,M)

if (toCalculate): 
    lam, vec, mesh, geometry = solve(geometry, thickness)
    results = s.convert_to_results(lam, vec, mesh, geometry)
    save_results(filename, results)
#    
#    save_mesh(meshfile, mesh)
#    
#    save_geometry(geometryfile, geometry)
    
else:
    results = load_results(filename)

    
    results_index = 0
    
#    plot.plot_mesh(results[results_index].mesh, width, thickness)
    
#    plot.plot_deformed_mesh(results[results_index], width, thickness)
    
#    plot.plot_init_and_deformed_geometry(results[results_index], 0, width, -thickness / 2, thickness / 2, 0)
    
#    plot.plot_init_geometry(results[results_index].geometry, 0, width, -thickness / 2, thickness / 2, 0)
    
    for i in range(6):
        plot.plot_strain_2(results[results_index], N, M, 0, width, -thickness / 2, thickness / 2, 0, i)
#        plot.plot_strain(results[results_index], 0, width, -thickness / 2, thickness / 2, 0, i)
    
    
    to_print = 20
    if (len(results) < to_print):
        to_print = len(results)
    
    for i in range(to_print):
        print(results[i].freq)
