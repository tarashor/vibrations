import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.general2D.solverlinear as s
import fem.general2D.result2D as r
import fem.mesh as me
import plot

from fem.general2D.matrices2D import stiffness_matrix, mass_matrix


def solve(geometry, layers, N, M):
    print(layers)
    
#    model = m.Model(geometry, layers, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
    model = m.Model(geometry, layers, m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)
    mesh = me.Mesh.generate2D(geometry.width, layers, N, M, model.boundary_conditions)
#    print(mesh.nodes)
    
    lam, vec = s.solve(model, mesh, stiffness_matrix, mass_matrix)
    
    return lam, vec, mesh, geometry



material1 = mat.IsotropicMaterial.steel()
material2 = mat.IsotropicMaterial.rubber()

width = 1
curvature = 0
thickness = 0.1

rubber_h_coef = 0.91
rh = rubber_h_coef*thickness
sh = (1-rubber_h_coef)*thickness/2

layers = set()
if (rubber_h_coef != 0 and rubber_h_coef != 1):
    rh = rubber_h_coef*thickness
    sh = (1-rubber_h_coef)*thickness/2
    
    h = sh+rh+sh
    top = h/2
    layer1 = m.Layer(top - sh, top, material1, 0)
    layer2 = m.Layer(top - sh-rh, top-sh, material2, 1)
    layer3 = m.Layer(top - sh-rh - sh, top - sh-rh, material1, 2)
    layers.add(layer1)
    layers.add(layer2)
    layers.add(layer3)
elif (rubber_h_coef == 0):
    sh = thickness
    
    top = thickness/2
    layer1 = m.Layer(top - sh, top, material1, 0)
    
    layers.add(layer1)
    
elif (rubber_h_coef == 1):
    rh = thickness
    
    top = thickness/2
    layer1 = m.Layer(top - rh, top, material2, 0)
    
    layers.add(layer1)

corrugation_amplitude = 0
corrugation_frequency = 0

geometry = g.General(width, curvature, corrugation_amplitude, corrugation_frequency)

N = 50
M = 2



lam, vec, mesh, geometry = solve(geometry, layers, N, M)
results = r.Result.convert_to_results(lam, vec, mesh, geometry)

results_index = 0

result = results[results_index]


plot.plot_init_and_deformed_geometry_in_cartesian(results[results_index], 0, width, -thickness / 2, thickness / 2, 0, geometry.to_cartesian_coordinates)

#plot.plot_deformed_mesh(results[results_index], width, thickness)

to_print = 20
if (len(results) < to_print):
    to_print = len(results)

for i in range(to_print):
    print(results[i].freqHz())

#tN = 1000
#T = 0.03
#
#
#print(result.freqHz())
#
#
#x = []
#y = []
#
#deltat = T / tN
#
#for t in range(tN):
#    time = deltat*t
#    u = result.get_displacement_and_deriv(width/2, 0, 0, time)
#    u1 = u[0]
#    u2 = u[4]
#    u3 = u[8]
#    x.append(t)
#    y.append(u3)
#    
#
#plot.plot_vibrations(x,y)
