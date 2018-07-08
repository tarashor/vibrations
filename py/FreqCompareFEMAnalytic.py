import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.mesh as me

import fem.shells1D.secondorder.shellsolver as s
import fem.shells1D.secondorder.result1D as r

#import fem.general2D.solverlinear as s
#import fem.general2D.result2D as r
#
#import fem.general2D.solver_nonlinear as s_nl
#import fem.general2D.result2Dnonlinear as r_nl

import plot

import numpy as np

#from fem.general2D.matrices2D import stiffness_matrix, mass_matrix, stiffness_matrix_nl_1, stiffness_matrix_nl_2
from fem.shells1D.secondorder.matrices1D import stiffness_matrix, mass_matrix, stiffness_matrix_nl_1, stiffness_matrix_nl_2




def solveLinear(geometry, thickness, material, N, M, u, bc):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, bc)
#    model = m.Model(geometry, layers, m.Model.FIXED_LEFT_RIGHT_EDGE)
    mesh = me.Mesh.generate1D(geometry.width, layers, N, model.boundary_conditions)
    
    lam, vec = s.solve(model, mesh, stiffness_matrix, mass_matrix)
    
    results_index = 0
    results = r.Result.convert_to_results(lam, vec, mesh, geometry, thickness)
    
    return results

#def solveLinear(geometry, thickness, material, N, M, u_max, bc):
#    layers = m.Layer.generate_layers(thickness, [material])
#    model = m.Model(geometry, layers, bc)
##    model = m.Model(geometry, layers, m.Model.FIXED_LEFT_RIGHT_EDGE)
#    mesh = me.Mesh.generate2D(geometry.width, layers, N, M, model.boundary_conditions)
#    
#    lam, vec = s.solve(model, mesh, stiffness_matrix, mass_matrix, u_max)
#    
##    results_index = 0
#    results = r.Result.convert_to_results(lam, vec, mesh, geometry)
#    
#    return results


def wAnalyticalLin(geometry, thickness, material, N, M, u_max):
    
    ko = 14/15
    G = material.C[4,4]
    LAM = ko * thickness * G
    
    c2 = np.sqrt(LAM/(rho*thickness))
    
    lam = np.pi / geometry.width
    
    alpha = (1+v)*v*v/(1-v-2*v*v)
    
    B = E*thickness/(1-v*v)*(1+alpha)
    
    D = (thickness*thickness/12)*B
    
    k12 = LAM / D 
    
    w = c2*lam*lam/(np.sqrt(k12+lam*lam))
    
    return w

def wAnalyticalLin2(geometry, thickness, material, bc):
    
    c = np.sqrt(E/rho)
    
    w = c*thickness*np.pi*np.pi/(np.sqrt(1-v*v)*geometry.width*geometry.width)
    
    if (bc == m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS):
        w = w*1/np.sqrt(12)
        print('points')
    else:
        w = w*2/3
        print('edges')
    
    return w


def wAnalyticalLin3(geometry, thickness, material, results_count):
    results = []
    
    h = thickness
    
    G = material.C[4,4]
    
    D = material.C[0,0]*h*h*h/12
    
    for i in range(1, results_count):
        a = i*np.pi/geometry.width
        
        res1 = G*a*(1+a)/(rho*h)
        res2 = (D*a*a+G*(1+a))/(rho*h)
        
        results.append(np.sqrt(res1)/(2*np.pi))
        results.append(np.sqrt(res2)/(2*np.pi))
    
    
    return results





E = 40*(10**9)
#E = 40000
v = 0.3
rho = 8000

kE3 = 100000000
#kE3 = 1
#kG13 = 100000000
kG13 = 1

material = mat.OrthotropicMaterial.create_from_E_and_v_with_koef_E3(E,v,rho)



#material.C[0,2] *= kE3
material.C[2,2] *= kE3
material.C[4,4] *= kG13

print(material.C)

#material = mat.IsotropicMaterial(E,v,rho)

width = 1
curvature = 0
thickness = 0.1

corrugation_amplitude = 0
corrugation_frequency = 0

geometry = g.General(width, curvature, corrugation_amplitude, corrugation_frequency)

N = 100
M = 4

#bc = m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS
bc = m.Model.FIXED_LEFT_RIGHT_EDGE

norm_koef = 2
u_max = norm_koef*thickness

results = solveLinear(geometry, thickness, material, N, M, u_max, bc)

wa = wAnalyticalLin(geometry, thickness, material, N, M, u_max)/(2*np.pi)
wa2 = wAnalyticalLin2(geometry, thickness, material, bc)/(2*np.pi)

results_count = 6

wa3 = wAnalyticalLin3(geometry, thickness, material, results_count)

for i in range(results_count):
    print('FEM {} = {}'.format(i, results[i].freqHz()))

print (wa3)

print('Analyt = {}'.format(wa))

print('Analyt2 = {}'.format(wa2))


