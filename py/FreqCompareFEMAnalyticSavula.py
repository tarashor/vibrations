import fem.geometry as g
import fem.model as m
import fem.material as mat
import fem.mesh as me

import fem.shells1D.KirchhoffLove.shellsolver as sKL
import fem.shells1D.KirchhoffLove.result1D as rKL
import fem.shells1D.KirchhoffLove.matrices1D as matKL

import fem.shells1D.firstorder.shellsolver as s1D1O
import fem.shells1D.firstorder.result1D as r1D1O
import fem.shells1D.firstorder.matrices1D as mat1D1O

import fem.shells1D.secondorder.shellsolver as s1D2O
import fem.shells1D.secondorder.result1D as r1D2O
import fem.shells1D.secondorder.matrices1D as mat1D2O

import fem.general2D.solverlinear as s2D
import fem.general2D.result2D as r2D
import fem.general2D.matrices2D as mat2D
#
#import fem.general2D.solver_nonlinear as s_nl
#import fem.general2D.result2Dnonlinear as r_nl

import plot

import matplotlib.pyplot as plt

import numpy as np


def solveLinearShellKirchhoffLove(geometry, thickness, material, N, bc):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, bc)
    mesh = me.Mesh.generate1D(geometry.width, layers, N, model.boundary_conditions)
    
    lam, vec = sKL.solve(model, mesh, matKL.stiffness_matrix, matKL.mass_matrix)
    
#    print(lam[1:20])
    
    results = rKL.Result.convert_to_results(lam, vec, mesh, geometry)
    
    return results

def solveLinearShell1Order(geometry, thickness, material, N, bc):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, bc)
    mesh = me.Mesh.generate1D(geometry.width, layers, N, model.boundary_conditions)
    
    lam, vec = s1D1O.solve(model, mesh, mat1D1O.stiffness_matrix, mat1D1O.mass_matrix)
#    print(lam)
    results = r1D1O.Result.convert_to_results(lam, vec, mesh, geometry)
    
    return results


def solveLinearShell2Order(geometry, thickness, material, N, bc):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, bc)
    mesh = me.Mesh.generate1D(geometry.width, layers, N, model.boundary_conditions)
    
    lam, vec = s1D2O.solve(model, mesh, mat1D2O.stiffness_matrix, mat1D2O.mass_matrix)
    
    results = r1D2O.Result.convert_to_results(lam, vec, mesh, geometry, thickness)
    
    return results

def solveLinear2D(geometry, thickness, material, N, M, bc):
    layers = m.Layer.generate_layers(thickness, [material])
    model = m.Model(geometry, layers, bc)
    mesh = me.Mesh.generate2D(geometry.width, layers, N, M, model.boundary_conditions)
    
    lam, vec = s2D.solve(model, mesh, mat2D.stiffness_matrix, mat2D.mass_matrix, u_max)
    
    results = r2D.Result.convert_to_results(lam, vec, mesh, geometry)
    
    return results


def wAnalyticalLin(geometry, thickness, material, N, M, u_max):
    
    ko = 5/6
    G = material.C[4,4]
    LAM = ko * thickness * G
    
    c2 = np.sqrt(LAM/(rho*thickness))
    
    lam = np.pi / geometry.width
    
    alpha = (1+v)*v*v/(1-v-2*v*v)
    
    B = E*thickness/(1-v*v)*(1+alpha)
    
    D = (thickness*thickness/12)*B
    
    k12 = LAM / D 
    
    w = 4/np.sqrt(3)*c2*lam*lam/(np.sqrt(k12+4*lam*lam))
    
    u = 2*lam*np.sqrt(B/((thickness*rho)))
    
    return w, u

def wAnalyticalLin2(geometry, thickness, material, bc):
    Dh = material.C[0,0]
    c = np.sqrt(Dh/rho)
    
    w = c*thickness*np.pi*np.pi/(geometry.width*geometry.width)
    
    if (bc == m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS):
        w = w*1/np.sqrt(12)
        print('points')
    else:
        w = w*2/3
        print('edges')
    
    return w


def wAnalyticalLin3(geometry, thickness, material, rGW, i):
    
    h = thickness
    
    G = material.C[4,4]*h
    
    print('G = {}'.format(G))
    
    h3 = h*h*h
    
    D = material.C[0,0]*h3/12
    print('D = {}'.format(D))
    
    lamb = np.pi/geometry.width
    print('Diff = {}'.format( (5*G*np.pi)/(geometry.width*(5*G+24*lamb*lamb*D))))
    
    a = i*np.pi/geometry.width
    
    res1 = G*a*(rGW+a)
    
    res2 = (D*a*a+G*(1+a/rGW))
        
    return np.sqrt(res1)/(2*np.pi), np.sqrt(res2)/(2*np.pi)
    



E = 40*(10**9)
#E = 1
v = 0.3
rho = 8000

#kE3 = 100000000
kE3 = 1
#kG13 = 100000000
kG13 = 1

material = mat.OrthotropicMaterial.create_from_E_and_v_with_koef_E3(E,v,rho, kE3)


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

M = 1

#bc = m.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS
bc = m.Model.FIXED_LEFT_RIGHT_EDGE

norm_koef = 2
u_max = norm_koef*thickness

results2D = solveLinear2D(geometry, thickness, material, N, M, bc)
results1D2O = solveLinearShell2Order(geometry, thickness, material, N, bc)
results1D1O = solveLinearShell1Order(geometry, thickness, material, N, bc)
resultsKL = solveLinearShellKirchhoffLove(geometry, thickness, material, N, bc)

wa, ua = wAnalyticalLin(geometry, thickness, material, N, M, u_max)
wa2 = wAnalyticalLin2(geometry, thickness, material, bc)/(2*np.pi)

results_count = len(resultsKL)

#wa3 = wAnalyticalLin3(geometry, thickness, material, results_count)

i = 0
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!CHANGE Mas matrix to 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print('KL {} = {}'.format(i, results1D1O[i].freqHz()))

Gni = np.argmax(np.absolute(resultsKL[i].g))
Wni = np.argmax(np.absolute(resultsKL[i].w))

x1s = sorted([n.x1 for n in resultsKL[i].mesh.nodes])

Gmax = resultsKL[i].g[Gni]/np.cos((i+1)*np.pi*x1s[Gni]/geometry.width)

Wmax = resultsKL[i].w[Wni]/np.sin((i+1)*np.pi*x1s[Wni]/geometry.width)

print(Gmax)
print(Wmax)
    
GW = Gmax/Wmax


x = sorted([n.index for n in results1D1O[i].mesh.nodes])
plt.figure()
plt.plot(x, results1D1O[i].u, 'y')
plt.plot(x, results1D1O[i].g, 'r')
plt.plot(x, results1D1O[i].w, 'g')
plt.grid()
plt.show()

res1, res2 = wAnalyticalLin3(geometry, thickness, material, GW, i+1)

print('res1 = {}'.format(res1))
print('res2 = {}'.format(res2))

print('anal wa, ua = {}'.format((wa/(2*np.pi), ua/(2*np.pi))))
#
#i = results_count//2
#print('KL {} = {}'.format(i, resultsKL[i].freqHz()))

#x = sorted([n.index for n in resultsKL[i].mesh.nodes])
#plt.figure()
#plt.plot(x, resultsKL[i].g, 'r')
#plt.plot(x, resultsKL[i].w, 'g')

#for i in range(results_count):
#    print('2D {} = {}'.format(i, results2D[i].freqHz()))
#    print('1D2O {} = {}'.format(i, results1D2O[i].freqHz()))
#    print('1D1O {} = {}'.format(i, results1D1O[i].freqHz()))
    
#    
    
#    print(resultsKL[i].g)
#    print(resultsKL[i].w)
    
#    print("i = {}".format(i))
#    
#    x = sorted([n.index for n in resultsKL[i].mesh.nodes])
#    
#    plt.figure()
#    
#    plt.plot(x, resultsKL[i].g, 'r')
#    
#    plt.plot(x, resultsKL[i].w, 'g')
#
#    
#    plt.grid()
#    plt.show()
    
    
#    Gni = np.argmax(np.absolute(resultsKL[i].g))
#    Wni = np.argmax(np.absolute(resultsKL[i].w))
#    
#    GW = resultsKL[i].g[Gni]/resultsKL[i].w[Wni]
#    print(GW)
#    
#    res1, res2 = wAnalyticalLin3(geometry, thickness, material, GW, i+1)
#    
#    print('KL {} = {}'.format(i, resultsKL[i].freqHz()))
#    print('res1 {} = {}'.format(i, res1))
#    print('res2 {} = {}'.format(i, res2))

#print (sorted(wa3))

#print('Analyt = {}'.format(wa))

#print('Analyt Volmir = {}'.format(wa2))


