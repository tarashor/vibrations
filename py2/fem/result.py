import numpy as np
from . import finiteelements as fe
from . import matrices
from math import cos


class Result:
    def __init__(self, freq, u1, u2, u3, mesh, geometry):
        self.freq = freq
        self.u1 = u1
        self.u2 = u2
        self.u3 = u3
        self.mesh = mesh
        self.geometry = geometry

    def get_displacement_and_deriv(self, x1, x2, x3, time):
        element = self.mesh.get_element(x1, x3)

        if (element is None):
            print ("x1 = {}, x2 = {}".format(x1, x3))

        u_nodes = np.zeros((8))

        u_nodes[0] = self.u1[element.top_left_index]
        u_nodes[1] = self.u1[element.top_right_index]
        u_nodes[2] = self.u1[element.bottom_right_index]
        u_nodes[3] = self.u1[element.bottom_left_index]

        u_nodes[4] = self.u3[element.top_left_index]
        u_nodes[5] = self.u3[element.top_right_index]
        u_nodes[6] = self.u3[element.bottom_right_index]
        u_nodes[7] = self.u3[element.bottom_left_index]

        h_e = matrices.element_aprox_functions(element, x1, x2, x3)

        return h_e.dot(u_nodes) * self.fi(time)

    def get_strain(self, x1, x2, x3, time):

        B = matrices.deriv_to_grad(self.geometry, x1, x2, x3)

        u = self.get_displacement_and_deriv(x1, x2, x3, time)

        grad_u = B.dot(u)

        E = matrices.grad_to_strain()
#        E_NL = grad_to_strain_nonlinear_matrix(alpha1, alpha2, geometry, grad_u)
        return E.dot(grad_u)

    def fi(self, time):
        return cos(self.freq * time)
