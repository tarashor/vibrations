class Model:

	SHARNIR = 1
	JORSTKO_LEFT = 2
	JORSTKO_LEFT_RIGHT = 3

	def __init__(self, geometry, layers, boundary_conditions):
		self.geometry = geometry
		self.layers = layers
		self.boundary_conditions = boundary_conditions

class Geometry:
	def __init__(self, width, curvature, wave_amp, wave_frequency):
		self.width = width
		self.curvature = curvature
		self.wave_amp = wave_amp
		self.wave_frequency = wave_frequency

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

	def height(self):
		return self.top - self.bottom

