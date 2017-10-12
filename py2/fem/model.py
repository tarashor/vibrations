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
    def __init__(self, width, curvature, corrugation_amplitude, corrugation_frequency):
        self.width = width
        self.curvature = curvature

    def __get_metric_tensor_components(self, alpha1, alpha2):
        q = 1 + self.curvature * alpha2
        a = self.curvature * (alpha1 - self.width/2)
        return q, a

    def get_g_11(self, alpha1, alpha2):
        q, a = self.__get_metric_tensor_components(alpha1, alpha2)
        return 1 / (q * q)

    def get_G(self, alpha1, alpha2):
        q, a = self.__get_metric_tensor_components(alpha1, alpha2)

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
