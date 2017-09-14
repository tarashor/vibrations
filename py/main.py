
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

class MeshNode(object):
	"""docstring for MeshNode"""
	def __init__(self, x, y):
		self.x = x
		self.y = y
		

class MeshElement(object):
	"""docstring for MeshElement"""
	def __init__(self, left_x, right_x, top_y, bottom_y, top_left_index, top_right_index, bottom_right_index, bottom_left_index):
		self.left_x = left_x
		self.top_y = top_y
		self.right_x = right_x
		self.bottom_y = bottom_y

		self.top_left_index = top_left_index
		self.top_right_index = top_right_index
		self.bottom_right_index = bottom_right_index
		self.bottom_left_index = bottom_left_index
		
	def height(self):
		return self.top_y - self.bottom_y

	def width(self):
		return self.right_x - self.left_x

	def __repr__(self):
		#return 'Element(x={}, y={})'.format(self.left_x, self.top_y)
		return 'Element(0="{}", 1="{}", 2="{}", 3="{}")'.format(self.top_left_index, self.top_right_index, self.bottom_right_index, self.bottom_left_index)

class Mesh(object):
	"""docstring for Mesh"""
	def __init__(self, elements):
		self.elements = elements

	@staticmethod
	def generate(width, layers, elements_width, elements_height_per_layer):
		d_x = width / elements_width
		# d_y = {}
		# for layer in layers:
		# 	d_y[layer] = layer.height() / elements_height_per_layer
		elements = []

		for layer in layers:
			d_y = layer.height() / elements_height_per_layer
			y = layer.top
			for i in range(elements_height_per_layer):
				x = 0
				for j in range(elements_width):
					top_left_index = i*(N+1)+j
					top_right_index =top_left_index + 1
					bottom_left_index = top_left_index + N + 1
					botton_right_index = bottom_left_index + 1
					element = MeshElement(x, x+d_x, y, y-d_y, top_left_index, top_right_index, botton_right_index, bottom_left_index)
					elements.append(element)
					x += d_x
				y -= d_y

		return Mesh(elements)

# def solve(model, mesh):
# 	pass


l = 1
K = 0.8
h = 0.1

gA=0.03
gV=20

N = 10
M = 2

geometry = Geometry(l,K, gA, gV)
layer = Layer(-h/2, h/2, Material.steel())
layers = [layer]
print(layers)
model = Model(geometry, layers, Model.SHARNIR)

mesh = Mesh.generate(model.geometry.width, layers, N, M)
print(mesh.elements)
print(len(mesh.elements))


# [vec lam] = solve(geom, layerModel, N, M, staticIndecies);

# ind = 1;

# printf ("Minimum frequancy = %f\n", sqrt(lam(ind)/rho));
# resVector = vec(:, ind);

# printf ("Norm of resVector = %f\n", sqrt(resVector'*resVector));

# midPaneResult=zeros(N+1,1);
# for p=1:N+1
#   i_new = (M+1)*(N+1)+(M/2)*(N + 1) + p;
#   temp=resVector(i_new);
#   midPaneResult(p) = temp;
# end

# w=sqrt(lam(ind)/rho);
# y = (cos(w*0)+sin(w*0)).*midPaneResult;
# x=0:l/N:l;
# h = plot(x, y);
# axis([0 l -1 1]);
# for t = 0.001:0.0001:5
#   y = (cos(w*t)+sin(w*t)).*midPaneResult;
#   set(h, 'YData', y);
#   pause(0.1);
# end