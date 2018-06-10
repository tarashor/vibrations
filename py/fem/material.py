#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 20:31:07 2018

@author: tarasgoriachko
"""

import numpy as np

class OrthotropicMaterial:
    def __init__(self, C11, C12, C13, C22, C23, C33, C44, C55, C66, rho):
        self.C = self.__create_matrix_C(C11, C12, C13, C22, C23, C33, C44, C55, C66)
        self.rho = rho
        
    def __create_matrix_C(self, C11, C12, C13, C22, C23, C33, C44, C55, C66):
        C = np.zeros((6, 6))

        C[0, 0] = C11
        C[1, 1] = C22
        C[2, 2] = C33
        C[0, 1] = C[1, 0] = C12
        C[0, 2] = C[2, 0] = C13
        C[1, 2] = C23
        C[2, 1] = C23

        C[3, 3] = C44
        C[4, 4] = C55
        C[5, 5] = C66

        return C


    def matrix_C(self, geometry, x1, x2, x3):
        return self.C
    
#    def create_from_E_and_v(E, v, G):
#        D = (1- v[1,2]*v[2,1]-v[2,3]*v[3,2]-v[1,3]*v[3,1]-2*v[1,3]*v[2,1]*v[3,2])/(E[0]*E[1]*E[2])
#        C11=(1-v[1,2]*v[2,1])/(D*E[0]*E[1])
#        C12=(1-v[1,2]*v[2,1])/(D*E[0]*E[1])
#        C13=(1-v[1,2]*v[2,1])/(D*E[0]*E[1])
#        C22=(1-v[1,2]*v[2,1])/(D*E[0]*E[1])
#        C23=(1-v[1,2]*v[2,1])/(D*E[0]*E[1])
#        C33=(1-v[1,2]*v[2,1])/(D*E[0]*E[1])
#        
#        return OrthotropicMaterial(C11, C12, C13, C22, C23, C33, G[0], G[1], G[2])
    

    
class IsotropicMaterial(OrthotropicMaterial):
    def __init__(self, E, v, rho):
        koef = E / ((1 + v) * (1 - 2 * v))
        C11=koef*(1 - v)
        C12=koef*v
        C44=koef*(1 - 2 * v) * 0.5
        super().__init__(C11, C12, C12, C11, C12, C11, C44, C44, C44, rho)
        
        
    def lam(self):
        return self.C[0,1]
    
    def mu(self):
        return self.C[3,3]
    
        
        

#    def tensor_C(self, geometry, x1, x2, x3):
#        C = np.zeros((6, 6))
#
#        v = self.v
#
#        C[0, 0] = (1 - v)
#        C[1, 1] = 1 - v
#        C[2, 2] = 1 - v
#        C[0, 1] = C[1, 0] = v
#        C[0, 2] = C[2, 0] = v
#        C[1, 2] = v
#        C[2, 1] = v
#
#        C[3, 3] = (1 - 2 * v) * 0.5
#        C[4, 4] = (1 - 2 * v) * 0.5
#        C[5, 5] = (1 - 2 * v) * 0.5
#
#        koef = self.E / ((1 + v) * (1 - 2 * v))
#
#        return koef * C

#    def tensor_C(self, geometry, x1, x2, x3):
#        N = 6
#
#        C = np.zeros((N, N))
#
#        lam = self.v * self.E / ((1 + self.v) * (1 - 2 * self.v))
#        mu = self.E / ((1 + self.v) * 2)
#
#        g = geometry.metric_tensor(x1, x2, x3)
#        
#        g = np.linalg.inv(g)
#
#        for i in range(N):
#            for j in range(N):
#                n, m = self.__get_index_conv(i)
#                k, l = self.__get_index_conv(j)
#                C[i, j] = mu * (g[n, k] * g[m, l] + g[n, l] * g[m, k]) + lam * g[n, m] * g[k, l]
#
#        return C
#
#    def __get_index_conv(self, index):
#        i = 0
#        j = 0
#        if (index == 0):
#            i = 0
#            j = 0
#        elif (index == 1):
#            i = 1
#            j = 1
#        elif (index == 2):
#            i = 2
#            j = 2
#        elif (index == 3):
#            i = 0
#            j = 1
#        elif (index == 4):
#            i = 0
#            j = 2
#        elif (index == 5):
#            i = 1
#            j = 2
#
#        return i, j

    @staticmethod
    def steel():
        return IsotropicMaterial(210000000000, 0.3, 8000)