import numpy as np
from . import matrices2D as matrices
from math import cos


class ResultNL:
    def __init__(self):
        pass
    
    def __init__(self, freq, u1, u2, u3,u11, u12, u13,u21, u22, u23,u31, u32, u33, mesh, geometry):
        self.freq = freq
        self.u1 = u1
        self.u2 = u2
        self.u3 = u3
        self.u11 = u11
        self.u12 = u12
        self.u13 = u13
        self.u21 = u21
        self.u22 = u22
        self.u23 = u23
        self.u31 = u31
        self.u32 = u32
        self.u33 = u33
        self.mesh = mesh
        self.geometry = geometry
        
    def freqHz(self):
        return self.freq/(2*np.pi)


    def get_displacement_and_deriv(self, x1, x2, x3, time, isPrint = False):
        element = self.mesh.get_element(x1, x3)

        if (element is None or isPrint):
            print ("x1 = {}, x3 = {}".format(x1, x3))
            

        u_nodes = np.zeros((8))

        u_nodes[0] = self.u1[element.top_left_index]
        u_nodes[1] = self.u1[element.top_right_index]
        u_nodes[2] = self.u1[element.bottom_right_index]
        u_nodes[3] = self.u1[element.bottom_left_index]

        u_nodes[4] = self.u3[element.top_left_index]
        u_nodes[5] = self.u3[element.top_right_index]
        u_nodes[6] = self.u3[element.bottom_right_index]
        u_nodes[7] = self.u3[element.bottom_left_index]
        
        u_nodes1 = np.zeros((8))

        u_nodes1[0] = self.u11[element.top_left_index]
        u_nodes1[1] = self.u11[element.top_right_index]
        u_nodes1[2] = self.u11[element.bottom_right_index]
        u_nodes1[3] = self.u11[element.bottom_left_index]

        u_nodes1[4] = self.u13[element.top_left_index]
        u_nodes1[5] = self.u13[element.top_right_index]
        u_nodes1[6] = self.u13[element.bottom_right_index]
        u_nodes1[7] = self.u13[element.bottom_left_index]
        
        u_nodes2 = np.zeros((8))

        u_nodes2[0] = self.u21[element.top_left_index]
        u_nodes2[1] = self.u21[element.top_right_index]
        u_nodes2[2] = self.u21[element.bottom_right_index]
        u_nodes2[3] = self.u21[element.bottom_left_index]

        u_nodes2[4] = self.u23[element.top_left_index]
        u_nodes2[5] = self.u23[element.top_right_index]
        u_nodes2[6] = self.u23[element.bottom_right_index]
        u_nodes2[7] = self.u23[element.bottom_left_index]
        
        u_nodes3 = np.zeros((8))

        u_nodes3[0] = self.u31[element.top_left_index]
        u_nodes3[1] = self.u31[element.top_right_index]
        u_nodes3[2] = self.u31[element.bottom_right_index]
        u_nodes3[3] = self.u31[element.bottom_left_index]

        u_nodes3[4] = self.u33[element.top_left_index]
        u_nodes3[5] = self.u33[element.top_right_index]
        u_nodes3[6] = self.u33[element.bottom_right_index]
        u_nodes3[7] = self.u33[element.bottom_left_index]
        
        u=u_nodes-u_nodes1-u_nodes2-u_nodes3
        

        h_e = matrices.element_aprox_functions(element, x1, x2, x3)

        return h_e.dot(u) * self.fi(time)+h_e.dot(u_nodes1)+h_e.dot(u_nodes2) * self.fi(2*time)+h_e.dot(u_nodes3) * self.fi(3*time)
    

    def fi(self, time):
        return cos(self.freq * time)
    
    
    @staticmethod
    def convert_to_result(l, U, U1, U2, U3, mesh, geometry):
        
        freq = np.sqrt(l)
        u1 = U[0:mesh.nodes_count()]
        u3 = U[mesh.nodes_count():2 * mesh.nodes_count()]
        u2 = np.zeros((mesh.nodes_count()))
        
        u11 = U1[0:mesh.nodes_count()]
        u13 = U1[mesh.nodes_count():2 * mesh.nodes_count()]
        u12 = np.zeros((mesh.nodes_count()))
        
        u21 = U2[0:mesh.nodes_count()]
        u23 = U2[mesh.nodes_count():2 * mesh.nodes_count()]
        u22 = np.zeros((mesh.nodes_count()))
        
        u31 = U3[0:mesh.nodes_count()]
        u33 = U3[mesh.nodes_count():2 * mesh.nodes_count()]
        u32 = np.zeros((mesh.nodes_count()))
        res = ResultNL(freq, u1, u2, u3,u11, u12, u13,u21, u22, u23,u31, u32, u33, mesh, geometry)
        
        return res
