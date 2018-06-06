import numpy as np
from . import finiteelements1D as fe
from . import matrices1D
from math import cos


class Result:
    def __init__(self):
        pass
    
    def __init__(self, freq, u, g, w, mesh, geometry):
        self.freq = freq
        self.u = u
        self.g = g
        self.w = w
        self.mesh = mesh
        self.geometry = geometry
        
    def rad_per_sec_to_Hz(self, rps):
        return rps/(2*np.pi)


    def get_displacement_and_deriv(self, x1, x2, x3, time):
        element = self.mesh.get_element(x1)

        if (element is None):
            print ("x1 = {}, x2 = {}".format(x1, x3))

        u_nodes = np.zeros((6))

        u_nodes[0] = self.u[element.start_index]
        u_nodes[1] = self.g[element.start_index]
        u_nodes[2] = self.w[element.start_index]

        u_nodes[3] = self.u[element.top_left_index]
        u_nodes[4] = self.g[element.top_right_index]
        u_nodes[5] = self.w[element.bottom_right_index]

        h_e = matrices1D.element_aprox_functions(element, x1, x2, x3)
        u3D = matrices1D.ugw_to_u1u3(x1, x2, x3)

        return u3D.dot(h_e.dot(u_nodes)) * self.fi(time)

    def get_strain(self, x1, x2, x3, time):

        B = matrices1D.deriv_to_grad(self.geometry, x1, x2, x3)

        u = self.get_displacement_and_deriv(x1, x2, x3, time)

        grad_u = B.dot(u)

        E = matrices.grad_to_strain()
#        E_NL = grad_to_strain_nonlinear_matrix(alpha1, alpha2, geometry, grad_u)
        return E.dot(grad_u)
    
    def get_strain_nl(self, x1, x2, x3, time):

        B = matrices1D.deriv_to_grad(self.geometry, x1, x2, x3)

        u = self.get_displacement_and_deriv(x1, x2, x3, time)

        grad_u = B.dot(u)

        E = matrices1D.grad_to_strain()
        
        E_NL = matrices1D.deformations_nl(self.geometry, grad_u, x1, x2, x3)
        
        return (E + E_NL).dot(grad_u)
    

    def fi(self, time):
        return cos(self.freq * time)
