import numpy as np
from . import finiteelements1D as fe
from . import matrices1D
from math import cos


class Result:
    def __init__(self):
        pass
    
    def __init__(self, freq, u10, u11, u12, u30,u31, u32, mesh, geometry, thickness):
        self.freq = freq
        self.u10 = u10
        self.u11 = u11
        self.u12 = u12
        self.u30 = u30
        self.u31 = u31
        self.u32 = u32
        self.thickness=thickness
        self.mesh = mesh
        self.geometry = geometry
        
    def rad_per_sec_to_Hz(self, rps):
        return rps/(2*np.pi)


    def get_displacement_and_deriv(self, x1, x2, x3, time):
        element = self.mesh.get_element(x1)

        if (element is None):
            print ("x1 = {}, x2 = {}".format(x1, x3))

        u_nodes = np.zeros((12))

        u_nodes[0] = self.u10[element.start_index]
        u_nodes[1] = self.u11[element.start_index]
        u_nodes[2] = self.u12[element.start_index]
        u_nodes[3] = self.u30[element.start_index]
        u_nodes[4] = self.u31[element.start_index]
        u_nodes[5] = self.u32[element.start_index]
        
        u_nodes[6] = self.u10[element.end_index]
        u_nodes[7] = self.u11[element.end_index]
        u_nodes[8] = self.u12[element.end_index]
        u_nodes[9] = self.u30[element.end_index]
        u_nodes[10] = self.u31[element.end_index]
        u_nodes[11] = self.u32[element.end_index]

        h_e = matrices1D.element_aprox_functions(element, x1, x2, x3)
        u3D = matrices1D.ugw_to_u1u3(x1, x2, x3, self.thickness)

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
