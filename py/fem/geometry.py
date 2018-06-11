import numpy as np


class General:
    def __init__(self, width, inverted_radius, corrugation_amplitude, corrugation_frequency):
        self.width = width
        self.inverted_radius = inverted_radius
        self.corrugation_amplitude = corrugation_amplitude
        self.corrugation_frequency = corrugation_frequency
        
    def __get_trig_args(self, x1, x2, x3):
        a1 = np.pi / 2 + (self.width / 2 - x1)* self.inverted_radius

        a2 = 2 * np.pi * x1 / self.width
        
        return a1, a2

#    def get_A_and_K(self, x1, x2, x3):
#        a1, a2 = self.__get_trig_args(x1, x2, x3)
#        
#        z = 2*self.corrugation_amplitude * self.corrugation_frequency *np.pi*np.sin(self.corrugation_frequency*a2)/self.width
#        w = 1 +  self.corrugation_amplitude * self.inverted_radius * np.cos(self.corrugation_frequency * a2)
#        
#        A = np.sqrt(w*w+z*z)
#        
#        q = w * self.inverted_radius +self.corrugation_amplitude*(2*np.pi*self.corrugation_frequency/self.width)*(2*np.pi*self.corrugation_frequency/self.width)*np.cos(self.corrugation_frequency*a2)
#        K=(q*w+2*z*z* self.inverted_radius)/(A*A*A)
#        
#        
#        return A,K
    
    
    def get_A_and_K(self, x1, x2, x3):
        a1, a2 = self.__get_trig_args(x1, x2, x3)
        
        z = self.corrugation_amplitude * self.corrugation_frequency * np.sin(self.corrugation_frequency*a1)*self.inverted_radius
        w = 1 +  self.corrugation_amplitude * self.inverted_radius * np.cos(self.corrugation_frequency * a1)
        
        A = np.sqrt(w*w+z*z)
        
        q = w * self.inverted_radius + self.corrugation_amplitude*self.corrugation_frequency*self.corrugation_frequency*self.inverted_radius*self.inverted_radius*np.cos(self.corrugation_frequency*a1)
        K=(q*w+2*z*z* self.inverted_radius)/(A*A*A)
        
        
        return A,K

    def metric_tensor(self, x1, x2, x3):
        A,K = self.get_A_and_K(x1, x2, x3)
        g = np.zeros((3, 3))
        g[0, 0] = A*A*(1+x3*K)*(1+x3*K)
        g[1, 1] = g[2, 2] = 1
        return g
    
    def metric_tensor_inv(self, x1, x2, x3):
        A,K = self.get_A_and_K(x1, x2, x3)
        g = np.zeros((3, 3))
        g[0, 0] = 1/(A*A*(1+x3*K)*(1+x3*K))
        g[1, 1] = g[2, 2] = 1
        return g



    def to_cartesian_coordinates(self, x1, x2, x3):
        x = x1
        y = x2
        z = x3
        
        if (self.inverted_radius > 0):
            a1, a2 = self.__get_trig_args(x1, x2, x3)
            x = (1 / self.inverted_radius + x3 + self.corrugation_amplitude * np.cos(self.corrugation_frequency * a1)) * np.cos(a1)
            z = (1 / self.inverted_radius + x3 + self.corrugation_amplitude * np.cos(self.corrugation_frequency * a1)) * np.sin(a1)

        return x, y, z

    
    
    def getJacobian(self, x1, x2, x3):
        A,K = self.get_A_and_K(x1, x2, x3)
        return A#*(1+x3*K)

    def __str__(self):
        return "L={}, 1/R={}, g_a={}, g_v={}".format(self.width, self.inverted_radius, self.corrugation_amplitude, self.corrugation_frequency)

class Plate:
    def __init__(self, width):
        self.width = width
        
    def metric_tensor(self, x1, x2, x3):
        g = np.zeros((3, 3))
        g[0, 0] = g[1, 1] = g[2, 2] = 1
        return g
    
    def metric_tensor_inv(self, x1, x2, x3):
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
    
    def getJacobian(self, x1, x2, x3):
        return 1

    def __str__(self):
        return r"L={}".format(self.width)


class CylindricalPlate(Plate):
    def __init__(self, width, curvature):
        super().__init__(width)
        self.curvature = curvature

    def __get_metric_tensor_components(self, x1, x2, x3):
        q = 1 + self.curvature * x3
        return q

    def metric_tensor(self, x1, x2, x3):
        q = self.__get_metric_tensor_components(x1, x2, x3)
        g = super().metric_tensor(x1, x2, x3)

        g[0, 0] = (q * q)
        return g
    
    def metric_tensor_inv(self, x1, x2, x3):
        q = self.__get_metric_tensor_components(x1, x2, x3)
        g = super().metric_tensor_inv(x1, x2, x3)

        g[0, 0] = 1/(q * q)
        return g

    def kristophel_symbols(self, x1, x2, x3):
        G = super().kristophel_symbols(x1, x2, x3)
        q = self.__get_metric_tensor_components(x1, x2, x3)

        G[0, 0, 2] = -self.curvature * q
        G[0, 2, 0] = self.curvature / q
        G[2, 0, 0] = self.curvature / q

        return G

    def to_cartesian_coordinates(self, x1, x2, x3):
        x = x1
        y = x2
        z = x3
        if (self.curvature > 0):
            ar = (np.pi + self.curvature * self.width) / 2 - x1 * self.curvature
            x = (1 / self.curvature + x3) * np.cos(ar)
            z = (1 / self.curvature + x3) * np.sin(ar)

        return x, y, z
    
    def R1(self, x1, x2, x3):
        x,y,z = super().R1(x1, x2, x3)
        q = self.__get_metric_tensor_components(x1, x2, x3)
        if (self.curvature > 0):
            ar = (np.pi + self.curvature * self.width) / 2 - x1 * self.curvature
            x = 1/q * np.sin(ar)
            z = -1/q * np.cos(ar)
        return x, y, z
    
    def R3(self, x1, x2, x3):
        x,y,z = super().R2(x1, x2, x3)
        
        if (self.curvature > 0):
            ar = (np.pi + self.curvature * self.width) / 2 - x1 * self.curvature
            x = np.cos(ar)
            z = np.sin(ar)
        return x, y, z
    
    def getJacobian(self, x1, x2, x3):
        return 1 + self.curvature * x3

    def __str__(self):
        r = super().__str__()
        return r+", K={}".format(self.curvature)


class CorrugatedCylindricalPlate(CylindricalPlate):

    def __init__(self, width, curvature, corrugation_amplitude, corrugation_frequency):
        super().__init__(width, curvature)
        self.corrugation_amplitude = corrugation_amplitude
        self.corrugation_frequency = corrugation_frequency

    def __get_metric_tensor_components(self, x1, x2, x3):
        q = 1 + self.curvature * x3
        a = (np.pi + self.curvature * self.width) / 2 - x1 * self.curvature
        w = q + self.corrugation_amplitude * self.curvature * np.cos(self.corrugation_frequency * a)
        z = self.corrugation_amplitude * self.corrugation_frequency * self.curvature * np.sin(self.corrugation_frequency * a)
        return q, a, w, z

    def metric_tensor(self, x1, x2, x3):
        q, a, w, z = self.__get_metric_tensor_components(x1, x2, x3)
        g = super().metric_tensor(x1, x2, x3)

        w1_2 = 1 / (w * w)

        g[0, 0] = w1_2
        g[0, 2] = -z * w1_2
        g[2, 0] = -z * w1_2
        g[2, 2] = (w * w + z * z) * w1_2
        return g

    def kristophel_symbols(self, x1, x2, x3):
        G = super().kristophel_symbols(x1, x2, x3)
        q, a, w, z = self.__get_metric_tensor_components(x1, x2, x3)

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
            x = (1 / self.curvature + x3 + self.corrugation_amplitude * np.cos(self.corrugation_frequency * ar)) * np.cos(ar)
            z = (1 / self.curvature + x3 + self.corrugation_amplitude * np.cos(self.corrugation_frequency * ar)) * np.sin(ar)

        return x, y, z

    
    def R1(self, x1, x2, x3):
        x = x1
        y = x2
        z = x3
        q, a, w, z = self.__get_metric_tensor_components(x1, x2, x3)
        if (self.curvature > 0):
            ar = (np.pi + self.curvature * self.width) / 2 - x1 * self.curvature
            x = w * np.sin(ar) + z * np.cos(ar)
            y = -w * np.cos(ar) + z * np.sin(ar)
        return x, y, z
    
    def getJacobian(self, x1, x2, x3):
        q, a, w, z = self.__get_metric_tensor_components(x1, x2, x3)
        return w

    def __str__(self):
        r = super().__str__()
        return r + ", g_a={}, g_v={}".format(self.curvature, self.corrugation_amplitude, self.corrugation_frequency)
