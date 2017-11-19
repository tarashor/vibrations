import numpy as np

class Model:
    FIXED_BOTTOM_LEFT_RIGHT_POINTS = 1
    FIXED_LEFT_EDGE = 2
    FIXED_LEFT_RIGHT_EDGE = 3

    def __init__(self, geometry, layers, boundary_conditions):
        self.geometry = geometry
        self.layers = layers
        self.boundary_conditions = boundary_conditions

class Material:
    def __init__(self, E, v, rho):
        self.E = E
        self.v = v
        self.rho = rho
        
    def tensor_C(self, geometry, x1, x2, x3):
        C = np.zeros((6, 6))
        
        v = self.v
    
        C[0, 0] = (1 - v)
        C[1, 1] = 1 - v
        C[2, 2] = 1 - v
        C[0, 1] = C[1, 0] = v
        C[0, 2] = C[2, 0] = v
        C[1, 2] = v
        C[2, 1] = v
    
        C[3, 3] = (1 - 2 * v) * 0.5
        C[4, 4] = (1 - 2 * v) * 0.5
        C[5, 5] = (1 - 2 * v) * 0.5
    
        koef = self.E / ((1 + v) * (1 - 2 * v))
    
        return koef * C
    
#    def tensor_C(self, geometry, x1, x2, x3):
#        N = 6
#        
#        C = np.zeros((N, N))
#        
#        lam = self.v * self.E/((1+self.v)*(1-2*self.v))
#        mu = self.E/((1+self.v)*2)
#        
#        g = geometry.metric_tensor(x1, x2, x3)
#        
#        for i in range(N):
#            for j in range(N):
#                n,m = self.__get_index_conv(i)
#                k,l = self.__get_index_conv(j)
#                C[i,j] = mu *(g[n,k]*g[m,l]+g[n,l]*g[m,k])+lam*g[n,m]*g[k,l]
#       
#        return C
    
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
