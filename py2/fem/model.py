import numpy as np


class Model:
    FIXED_BOTTOM_LEFT_RIGHT_POINTS = 1
    FIXED_LEFT_EDGE = 2
    FIXED_LEFT_RIGHT_EDGE = 3

    def __init__(self, geometry, layers, boundary_conditions):
        self.geometry = geometry
        self.layers = layers
        self.boundary_conditions = boundary_conditions


class Geometry:
    def __init__(self, width, curvature):
        self.width = width
        self.curvature = curvature

    def __get_metric_tensor_components(self, alpha1, alpha2):
        q = 1 + self.curvature * alpha2
        return q

    def get_metric_tensor(self, alpha1, alpha2):
        q = self.__get_metric_tensor_components(alpha1, alpha2)
        G = np.zeros((3, 3))
        G[0, 0] = 1/(q*q)
        G[1, 1] = 1
        G[2, 2] = 1

        return G

    def get_G(self, alpha1, alpha2):
        q = self.__get_metric_tensor_components(alpha1, alpha2)

        G111 = 0
        G211 = -self.curvature * q
        G112 = self.curvature / q
        G121 = G112

        return G111, G211, G112, G121


class Material:
    def __init__(self, E, v, rho):
        self.E = E
        self.v = v
        self.rho = rho
        
#    def tensor_C1(self, alpha1, alpha2, geometry):
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
    
    def tensor_C(self, alpha1, alpha2, geometry):
        N = 6
        
        C = np.zeros((N, N))
        
        lam = self.v * self.E/((1+self.v)*(1-2*self.v))
        mu = self.E/((1+self.v)*2)
        
        g = geometry.get_metric_tensor(alpha1, alpha2)
        
        for i in range(N):
            for j in range(N):
                n,m = self.__get_index_conv(i)
                k,l = self.__get_index_conv(j)
                C[i,j] = mu *(g[n,k]*g[m,l]+g[n,l]*g[m,k])+lam*g[n,m]*g[k,l]
       
        return C
    
    def __get_index_conv(self, index):
        i = 0
        j = 0
        if (index == 0):
            i = 0
            j = 0
        elif (index == 1):
            i = 1
            j = 1
        elif (index == 2):
            i = 2
            j = 2
        elif (index == 3):
            i = 0
            j = 1
        elif (index == 4):
            i = 0
            j = 2
        elif (index == 5):
            i = 1
            j = 2
    
        return i, j

    @staticmethod
    def steel():
        return Material(210000000000, 0.3, 8000)


class Layer:
    def __init__(self, bottom, top, material, index_from_top):
        self.bottom = bottom
        self.top = top
        self.material = material
        self.index_from_top = index_from_top

    def __repr__(self):
        return "Layer({:f};{:f})".format(self.bottom, self.top)

    def __hash__(self):
        return hash(self.top)

    def __eq__(self, other):
        if (isinstance(other, Layer)):
            return self.top == other.top and self.height() == other.height() and self.index_from_top == other.index_from_top
        else:
            return False

    def height(self):
        return self.top - self.bottom
