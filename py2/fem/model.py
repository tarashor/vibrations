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
