import numpy as np
from . import matrices1D
from math import cos


class Result:
    def __init__(self):
        pass
    
    def __init__(self, freq, g, w, mesh, geometry):
        self.freq = freq
        self.g = g
        self.w = w
        self.mesh = mesh
        self.geometry = geometry
        
    @staticmethod
    def convert_to_results(eigenvalues, eigenvectors, mesh, geometry):
    
        results = []
        for i in range(eigenvalues.size):
            freq = np.sqrt(eigenvalues[i])
            g = eigenvectors[:, i][range(0,2 * mesh.nodes_count(),2)]
            w = eigenvectors[:, i][range(1,2 * mesh.nodes_count(),2)]
            r = Result(freq, g, w, mesh, geometry)
            results.append(r)
    
        return results
        
    def freqHz(self):
        return self.freq/(2*np.pi)
    
    def __u_to_u1u3(self, x1, x2, x3):
        T=np.zeros((12,4))
        T[0,0]=x3
        T[1,1]=x3
        T[3,0]=1
        
        T[8,2]=1
        T[9,3]=1
        return T


    def get_displacement_and_deriv(self, x1, x2, x3, time):
        element = self.mesh.get_element(x1, x2)

        if (element is None):
            print ("x1 = {}, x2 = {}".format(x1, x3))

        u_nodes = np.zeros((4))

        u_nodes[0] = self.g[element.start_index]
        u_nodes[1] = self.w[element.start_index]

        u_nodes[2] = self.g[element.end_index]
        u_nodes[3] = self.w[element.end_index]

        h_e = matrices1D.element_aprox_functions(element, x1, x2, x3)
        u3D = self.__u_to_u1u3(x1, x2, x3)

        return u3D.dot(h_e.dot(u_nodes)) * self.fi(time)

#    def get_strain(self, x1, x2, x3, time):
#
#        B = matrices1D.deriv_to_grad(self.geometry, x1, x2, x3)
#
#        u = self.get_displacement_and_deriv(x1, x2, x3, time)
#
#        grad_u = B.dot(u)
#
#        E = matrices.grad_to_strain()
##        E_NL = grad_to_strain_nonlinear_matrix(alpha1, alpha2, geometry, grad_u)
#        return E.dot(grad_u)
#    
#    def get_strain_nl(self, x1, x2, x3, time):
#
#        B = matrices1D.deriv_to_grad(self.geometry, x1, x2, x3)
#
#        u = self.get_displacement_and_deriv(x1, x2, x3, time)
#
#        grad_u = B.dot(u)
#
#        E = matrices1D.grad_to_strain()
#        
#        E_NL = matrices1D.deformations_nl(self.geometry, grad_u, x1, x2, x3)
#        
#        return (E + E_NL).dot(grad_u)
    

    def fi(self, time):
        return cos(self.freq * time)
