# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 13:53:13 2018

@author: Taras
"""

def contraction(A,B):
    res = A[0,0]*B[0,0]
    for i in range(3):
        for j in range(3):
            if (i != 0 or j != 0):
                res += A[i,j]*B[j,i]
    return res

T = zeros(3)
for i in range(3):
    for j in range(i, 3):
        T[i,j] = Symbol(r't_{{{}{}}}'.format(i+1, j+1))
        T[j,i] = T[i,j]
        
S = zeros(3)
for i in range(3):
    for j in range(i, 3):
        S[i,j] = Symbol(r's_{{{}{}}}'.format(i+1, j+1))
        S[j,i] = S[i,j]

A = zeros(3)
A[0,0]=Symbol(r'a_{11}')
A[1,1]=Symbol(r'a_{22}')
A[2,2]=Symbol(r'a_{33}')

G = zeros(3)
G[0,0]=A[0,0]**2
G[1,1]=A[1,1]**2
G[2,2]=A[2,2]**2

A_inv = A**-1

T_new = A_inv*T*A_inv.T
S_new = A_inv*S*A_inv.T

e=contraction(G*T_new*G,S_new)

display(e)