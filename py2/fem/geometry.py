import numpy as np


class Geometry:
    def metric_tensor(self, x1, x2, x3):
        g = np.zeros((3, 3))
        g[0, 0] = g[1, 1] = g[2, 2] = 1
        return g

    def kristophel_symbols(self, x1, x2, x3):
        G = np.zeros((3, 3, 3))
        return G

    def to_cartesian_coordinates(self, x1, x2, x3):
        return x1, x2, x3
    
    def R1(self, x1, x2, x3):
        return 1, 0, 0
    
    def R2(self, x1, x2, x3):
        return 0, 1, 0
    
    def R3(self, x1, x2, x3):
        return 0, 0, 1
    

    def __str__(self):
        return ""


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

        g[0, 0] = (q * q)
        return g

    def kristophel_symbols(self, x1, x2, x3):
        G = super().kristophel_symbols(x1, x2, x3)
        q = self.__get_metric_tensor_components(x1, x2)

        G[1, 0, 0] = -self.curvature * q
        G[0, 0, 1] = self.curvature / q
        G[0, 1, 0] = self.curvature / q

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
    
    def R1(self, x1, x2, x3):
        x,y,z = super().R1(x1, x2, x3)
        q = self.__get_metric_tensor_components(x1, x2)
        if (self.curvature > 0):
            ar = (np.pi + self.curvature * self.width) / 2 - x1 * self.curvature
            x = 1/q * np.cos(ar)
            y = 1/q * np.sin(ar)
        return x, y, z
    
    def R2(self, x1, x2, x3):
        x,y,z = super().R2(x1, x2, x3)
        
        if (self.curvature > 0):
            ar = (np.pi + self.curvature * self.width) / 2 - x1 * self.curvature
            x = np.sin(ar)
            y = np.cos(ar)
        return x, y, z
    

    def __str__(self):
        return "K={}".format(self.curvature)


class CorrugatedCylindricalPlate(CylindricalPlate):

    def __init__(self, width, curvature, corrugation_amplitude, corrugation_frequency):
        super().__init__(width, curvature)
        self.corrugation_amplitude = corrugation_amplitude
        self.corrugation_frequency = corrugation_frequency

    def __get_metric_tensor_components(self, x1, x2):
        q = 1 + self.curvature * x2
        a = (np.pi + self.curvature * self.width) / 2 - x1 * self.curvature
        w = q + self.corrugation_amplitude * self.curvature * np.cos(self.corrugation_frequency * a)
        z = self.corrugation_amplitude * self.corrugation_frequency * self.curvature * np.sin(self.corrugation_frequency * a)
        return q, a, w, z

    def metric_tensor(self, x1, x2, x3):
        q, a, w, z = self.__get_metric_tensor_components(x1, x2)
        g = super().metric_tensor(x1, x2, x3)

        w1_2 = 1 / (w * w)

        g[0, 0] = w1_2
        g[0, 1] = -z * w1_2
        g[1, 0] = -z * w1_2
        g[1, 1] = (w * w + z * z) * w1_2
        return g

    def kristophel_symbols(self, x1, x2, x3):
        G = super().kristophel_symbols(x1, x2, x3)
        q, a, w, z = self.__get_metric_tensor_components(x1, x2)

        G[0, 0, 0] = 2 * z * self.curvature / w
        G[1, 0, 0] = -self.curvature * self.curvature * self.corrugation_amplitude * self.corrugation_frequency * self.corrugation_frequency * np.cos(self.corrugation_frequency * a) - w * self.curvature - 2 * z * z * self.curvature / w
        
        G[0, 0, 1] = self.curvature / w
        G[0, 1, 0] = self.curvature / w
        
        G[1, 0, 1] = -z*self.curvature / w
        G[1, 1, 0] = -z*self.curvature / w

        return G

    def to_cartesian_coordinates(self, x1, x2, x3):
        x = x1
        y = x2
        z = x3
        
        if (self.curvature > 0):
            ar = (np.pi + self.curvature * self.width) / 2 - x1 * self.curvature
            x = (1 / self.curvature + x2 + self.corrugation_amplitude * np.cos(self.corrugation_frequency * ar)) * np.cos(ar)
            y = (1 / self.curvature + x2 + self.corrugation_amplitude * np.cos(self.corrugation_frequency * ar)) * np.sin(ar)

        return x, y, z

    
    def R1(self, x1, x2, x3):
        x = x1
        y = x2
        z = x3
        q, a, w, z = self.__get_metric_tensor_components(x1, x2)
        if (self.curvature > 0):
            ar = (np.pi + self.curvature * self.width) / 2 - x1 * self.curvature
            x = w * np.sin(ar) + z * np.cos(ar)
            y = -w * np.cos(ar) + z * np.sin(ar)
        return x, y, z
    

    def __str__(self):
        return "K={}, g_A={}, g_v={}".format(self.curvature, self.corrugation_amplitude, self.corrugation_frequency)
