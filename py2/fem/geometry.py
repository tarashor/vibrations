import numpy as np

class Geometry:
    def metric_tensor(self, x1, x2, x3):
        g = np.zeros((3,3))
        g[0,0] = g[1,1] = g[2,2] = 1
        return g


    def kristophel_symbols(self, x1, x2, x3):
        G = np.zeros((3,3,3))
        return G
    
    def to_cartesian_coordinates(self, x1, x2, x3):
        return x1, x2, x3



class CylindricalPlate(Geometry):
    def __init__(self, width, curvature):
        self.width = width
        self.curvature = curvature

    def __get_metric_tensor_components(self, x1, x2):
        q = 1 + self.curvature * x2
        return q

    def metric_tensor(self, x1, x2, x3):
        q = self.__get_metric_tensor_components(x1, x2)
        g = super().metric_tensor(x1, x2, x3)
        g[0,0] = 1 / (q * q)
        return g

    def kristophel_symbols(self, x1, x2, x3):
        G = super().kristophel_symbols(x1, x2, x3)
        q = self.__get_metric_tensor_components(x1, x2)

        G[1,0,0] = -self.curvature * q
        G[0,0,1] = self.curvature / q
        G[0,1,0] = self.curvature / q

        return G
    
    def to_cartesian_coordinates(self, x1, x2, x3):
        x = x1
        y = x2
        z = x3
        if (self.curvature > 0):
            ar = (np.pi + self.curvature * self.width) / 2 - x1 * self.curvature
            x = (1 / self.curvature + x2) * np.cos(ar)
            y = (1 / self.curvature + x2) * np.sin(ar)
            
        return x, y, z


class CorrugatedCylindricalPlate(CylindricalPlate):

    def __init__(self, width, curvature, corrugation_amplitude, corrugation_frequency):
        super().__init__(width, curvature)
        self.corrugation_amplitude = corrugation_amplitude
        self.corrugation_frequency = corrugation_frequency
