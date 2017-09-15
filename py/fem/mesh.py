from . import model

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
		
		elements = []

		for layer in layers:
			d_y = layer.height() / elements_height_per_layer
			y = layer.top
			for i in range(elements_height_per_layer):
				x = 0
				for j in range(elements_width):
					top_left_index = i*(elements_width+1)+j
					top_right_index =top_left_index + 1
					bottom_left_index = top_left_index + elements_width + 1
					botton_right_index = bottom_left_index + 1
					element = MeshElement(x, x+d_x, y, y-d_y, top_left_index, top_right_index, botton_right_index, bottom_left_index)
					elements.append(element)
					x += d_x
				y -= d_y

		return Mesh(elements)