import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.mesh as me


import fem.shells1D.secondorder.shellsolver as s1D2O
import fem.shells1D.secondorder.result1D as r1D2O
import fem.shells1D.secondorder.matrices1D as mat1D2O


def solveLinear1D2O(geometry, layers, N, bc):
    model = m.Model(geometry, layers, bc)
#    model = m.Model(geometry, layers, m.Model.FIXED_LEFT_RIGHT_EDGE)
    mesh = me.Mesh.generate1D(geometry.width, layers, N, model.boundary_conditions)
    
    lam, vec = s1D2O.solve(model, mesh, mat1D2O.stiffness_matrix, mat1D2O.mass_matrix)
    
    results_index = 0
    results = r1D2O.Result.convert_to_results(lam, vec, mesh, geometry, thickness)
    
    return results[results_index]



width = 1

corrugation_amplitude = 0
corrugation_frequency = 0

E = 40*(10**9)
v = 0.25
rho = 8000

#kE3 = 100000000
kE1 = 40
#kG13 = 100000000
kG13 = 1

material = mat.OrthotropicMaterial.create_from_E_and_v_with_koef_E1(E,v,rho, kE1)

material.C[4,4] *= kG13


N = 50

bc = m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS
#bc = m.Model.FIXED_LEFT_RIGHT_EDGE

thicknesses = [0.01, 0.05, 0.1, 0.2]

curvatures_coef = range(1,200,20)

y_per_hr = {}
x_per_hr = {}

for thickness in thicknesses:
    layers = m.Layer.generate_layers(thickness, [material])
    
    print('========== H = {} ============='.format(thickness))
    
    for cc in curvatures_coef:
        
        curvature = cc*thickness
        geometry = g.General(width, curvature, corrugation_amplitude, corrugation_frequency)
        
        
        result2D = solveLinear1D2O(geometry, layers, N, bc)
        
        print('--------------------------------------')
        print(cc)
        print(curvature)
        print(result2D.freqHz())
        print('--------------------------------------')
    
    
    print('========================')
        
    
        
    
    

