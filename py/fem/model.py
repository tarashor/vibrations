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
        self.corrugation_amplitude = corrugation_amplitude
        self.corrugation_frequency = corrugation_frequency


class Material:
    def __init__(self, E, v, rho):
        self.E = E
        self.v = v
        self.rho = rho

    @staticmethod
    def steel():
        return Material(210000000000, 0.3, 8000)


class Layer:
    def __init__(self, bottom, top, material):
        self.bottom = bottom
        self.top = top
        self.material = material

    def __repr__(self):
        return "Layer({};{})".format(self.bottom, self.top)

    def __hash__(self):
        return hash(self.top)

    def __eq__(self, other):
        if (isinstance(other, Layer)):
            return self.top == other.top and self.height() == other.height()
        else:
            return False

    def height(self):
        return self.top - self.bottom
