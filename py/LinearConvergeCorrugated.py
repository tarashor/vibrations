import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.mesh as me


import fem.shells1D.secondorder.shellsolver as s1D2O
import fem.shells1D.secondorder.result1D as r1D2O
import fem.shells1D.secondorder.matrices1D as mat1D2O

import plot



def solveLinear1D2O(geometry, layers, N, bc):
    model = m.Model(geometry, layers, bc)
    mesh = me.Mesh.generate1D(geometry.width, layers, N, model.boundary_conditions)
    
    lam, vec = s1D2O.solve(model, mesh, mat1D2O.stiffness_matrix, mat1D2O.mass_matrix)
    
    results_index = 0
    results = r1D2O.Result.convert_to_results(lam, vec, mesh, geometry, thickness)
    
#    for i in range(6):
#        plot.plot_init_and_deformed_geometry_in_cartesian(results[i], 0, width, -thickness / 2, thickness / 2, 0, geometry.to_cartesian_coordinates)
    
    return results[results_index]



width = 2
curvature = 0.8
thickness = 0.05

corrugation_amplitude = 0.03

material = mat.IsotropicMaterial.steel()

layers = m.Layer.generate_layers(thickness, [material])

N = 200

#bc = m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS
bc = m.Model.FIXED_LEFT_RIGHT_EDGE

#corrugation_frequencies = range(2, 50, 2)
corrugation_frequencies = [2, 4, 6, 8, 15, 20, 50, 80, 100, 200, 300, 500]
#corrugation_frequencies = [10, 20]
#corrugation_frequencies = [50]

v = []


for corrugation_frequency in corrugation_frequencies:

    geometry = g.General(width, curvature, corrugation_amplitude, corrugation_frequency)
        
    result1D = solveLinear1D2O(geometry, layers, N, bc)
    
#    plot.plot_init_and_deformed_geometry_in_cartesian(result1D, 0, width, -thickness / 2, thickness / 2, 0, geometry.to_cartesian_coordinates)
    
    f = result1D.freqHz()
    
    v.append(f)
    
    print('gv = {}, freq = {}'.format(corrugation_frequency, f))
    
    
#plot.plot_freq_from_corrugated_freq(corrugation_frequencies, v)
    
        
    
    

