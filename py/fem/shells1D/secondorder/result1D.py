import numpy as np
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
        
    def freqHz(self):
        return self.freq/(2*np.pi)
    
    def __u_to_u1u3(self, x1, x2, x3, h):
        p0=0.5-x3/h
        p1=0.5+x3/h
        p2=1-(2*x3/h)*(2*x3/h)
        
        dp0 = -1/h
        dp1 = 1/h
        dp2 = -8*x3/(h*h)
        
        L=np.zeros((12,12))
        
        L[0,0]=p0
        L[0,2]=p1
        L[0,4]=p2
        
        L[1,1]=p0
        L[1,3]=p1
        L[1,5]=p2
        
        L[3,0]=dp0
        L[3,2]=dp1
        L[3,4]=dp2
        
        L[8,6]=p0
        L[8,8]=p1
        L[8,10]=p2
        
        L[9,7]=p0
        L[9,9]=p1
        L[9,11]=p2
        
        L[11,6]=dp0
        L[11,8]=dp1
        L[11,10]=dp2
    
        return L
    
    @staticmethod
    def convert_to_results(eigenvalues, eigenvectors, mesh, geometry, thickness):
        results = []
        for i in range(eigenvalues.size):
            freq = np.sqrt(eigenvalues[i])
            u10 = eigenvectors[:, i][range(0,6 * mesh.nodes_count(),6)]
            u11 = eigenvectors[:, i][range(1,6 * mesh.nodes_count(),6)]
            u12 = eigenvectors[:, i][range(2,6 * mesh.nodes_count(),6)]
            u30 = eigenvectors[:, i][range(3,6 * mesh.nodes_count(),6)]
            u31 = eigenvectors[:, i][range(4,6 * mesh.nodes_count(),6)]
            u32 = eigenvectors[:, i][range(5,6 * mesh.nodes_count(),6)]
            r = Result(freq, u10, u11, u12, u30, u31, u32, mesh, geometry, thickness)
            results.append(r)
    
        return results


    def get_displacement_and_deriv(self, x1, x2, x3, time):
        element = self.mesh.get_element(x1, x2)

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
        u3D = self.__u_to_u1u3(x1, x2, x3, self.thickness)

        return u3D.dot(h_e.dot(u_nodes)) * self.fi(time)

    

    def fi(self, time):
        return cos(self.freq * time)
