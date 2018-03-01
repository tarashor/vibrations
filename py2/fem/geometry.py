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
    
    def normal_to_middle_surface(self, x1, x2, x3):
        return 0, 1, 0

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

        g[0, 0] = 1 / (q * q)
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
    
    def normal_to_middle_surface(self, x1, x2, x3):
        n1 = 0
        n2 = 1
        n3 = 0
        if (self.curvature > 0):
            ar = self.curvature * (self.width / 2 - x1)
            
            dr1 = np.cos(ar)
            dr2 = np.sin(ar)
            
            n1 = -dr2/np.sqrt(dr1*dr1+dr2*dr2)
            n2 = dr1/np.sqrt(dr1*dr1+dr2*dr2)

        return n1, n2, n3


    def __str__(self):
        return "K={}m".format(self.curvature)


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
        G[0, 0, 1] = self.curvature / q
        G[0, 1, 0] = self.curvature / q

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
    
    def normal_to_middle_surface(self, x1, x2, x3):
        n1 = 0
        n2 = 1
        n3 = 0
        if (self.curvature > 0):
            ar = (np.pi + self.curvature * self.width) / 2 - x1 * self.curvature
            
#            a=self.curvature*self.corrugation_frequency*self.corrugation_amplitude*np.sin(self.corrugation_frequency*ar)
#            b=(1+self.curvature*self.corrugation_amplitude*np.cos(self.corrugation_frequency*ar))
#
#            n1=-a*np.sin(ar)+b*np.cos(ar)
#            n2=a*np.cos(ar)+b*np.sin(ar)
#            
#            length_n=np.sqrt(a*a+b*b)
#
#
#            n1=n1/length_n
#            n2=n2/length_n
            
            dr1 = self.curvature*self.corrugation_frequency*self.corrugation_amplitude*np.cos(ar)*np.sin(self.corrugation_frequency*ar) + self.curvature*(self.corrugation_amplitude*np.cos(self.corrugation_frequency*ar)+1/self.curvature)*np.sin(ar)
            
            dr2 = self.curvature*self.corrugation_frequency*self.corrugation_amplitude*np.sin(ar)*np.sin(self.corrugation_frequency*ar) - self.curvature*(self.corrugation_amplitude*np.cos(self.corrugation_frequency*ar)+1/self.curvature)*np.cos(ar)
            
            n1 = -dr2/np.sqrt(dr1*dr1+dr2*dr2)
            n2 = dr1/np.sqrt(dr1*dr1+dr2*dr2)                   

        return n1, n2, n3
        
#        ar = pi/2+(L-2*x1)/(2*R);
#        r(1,:)=(R + g_a.*cos(g_f.*ar)).*cos(ar);
#r(2,:)=(R + g_a.*cos(g_f.*ar)).*sin(ar);
#
#a=g_a*g_f*K.*sin(g_f.*ar);
#b=(1+g_a*K.*cos(g_f.*ar));
#
#n(1,:)=-a.*sin(ar)+b.*cos(ar);
#n(2,:)=a.*cos(ar)+b.*sin(ar);
#
#length_n=sqrt(a.*a+b.*b); 
#
#n=n./length_n;

    def __str__(self):
        return r"$K={}$m, $g_A={}$m, $g_v={}$".format(self.curvature, self.corrugation_amplitude, self.corrugation_frequency)
