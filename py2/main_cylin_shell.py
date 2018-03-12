# -*- coding: utf-8 -*-

import fem.geometry as g
import fem.model as m
import fem.solver as s
import fem.mesh as me
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


def solve(width, curvature, thickness, Gkoef):
    layers_count = 1
    mat = m.Material.steel()
    mat.Gkoef = Gkoef
    layers = generate_layers(thickness, layers_count, mat)
    mesh = me.Mesh.generate(width, layers, N, M, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
    geometry = g.CylindricalPlate(width, curvature)
    # geometry = g.Geometry()
    model = m.Model(geometry, layers, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
    return s.solve(model, mesh)


# r=2
# width = r*2*3.14
# curvature = 1/r

width = 2
curvature = 0.8
thickness = 0.05


N = 50
M = 4

koefs = range(1,30)

w=[]

for Gkoef in koefs:
    results = solve(width, curvature, thickness, Gkoef)
    results_index = 0
    w.append(results[results_index].freq)
    print('{} - {}'.format(Gkoef, results[results_index].freq))
#    filename = "g_v {}".format(cf)
#    filename_def = "deform g_v {}".format(cf)
#    print({})
    #plot.plot_init_and_deformed_geometry(results[results_index], 0, width, -thickness / 2, thickness / 2, 0, filename_def)
#    plot.plot_init_geometry(results[results_index].geometry, 0, width, -thickness / 2, thickness / 2, filename)
#    plot.plot_normals(results[results_index].geometry, 0, width, -thickness / 2, thickness / 2)
# plot.plot_strain(results[results_index], 0, width, -thickness / 2, thickness / 2, 0)
    
    
#g = [2,4,6,8,10, 20, 50, 80, 100, 200, 300, 500]
#w = [799,775,692,1025,1056, 3658, 6383, 6936, 7709, 5559, 4914, 3195]
    
#g = [2,4,6,8,10, 20, 50, 80, 100]
#w = [799,775,692,1025,1056, 3658, 6383, 6936, 7709]
#
#plot.plot_freq_from_corrugated_freq(g, w, N, M)
    
plot.plot_2D_depen(koefs, w, r'$k$', r'$w_{min}$, Hz')


