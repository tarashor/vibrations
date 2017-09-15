import numpy as np
from numpy import linalg as la

def solve(model, mesh):
	s = stiffness_matrix(model, mesh)
	m = mass_matrix(model, mesh)

	s = mesh.apply_boundary_conditions(s)
	m = mesh.apply_boundary_conditions(m)
	lam, vec = la.eig(s,m)
	return lam

def stiffness_matrix(model, mesh):
	return np.zeros((2,2))

def mass_matrix(model, mesh):
	return np.zeros((2,2))

# b=np.array([[1,2], [3,4]])
# print(a.dot(b)==b)