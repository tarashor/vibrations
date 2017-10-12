from . import model


class MeshNode(object):
    """docstring for MeshNode"""

    def __init__(self, x, y, index):
        self.x = x
        self.y = y
        self.index = index

    def __eq__(self, other):
        return self.index == other.index

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        return "Node('{} - ({:f};{:f})')".format(self.index, self.x, self.y)


class MeshElement(object):
    """docstring for MeshElement"""

    def __init__(self, top_left, top_right, bottom_right, bottom_left):
        self.top_left = top_left
        self.top_right = top_right
        self.bottom_right = bottom_right
        self.bottom_left = bottom_left

    def height(self):
        return self.top_right.y - self.bottom_left.y

    def width(self):
        return self.top_right.x - self.bottom_left.x

    @property
    def top_left_index(self):
        return self.top_left.index

    @property
    def top_right_index(self):
        return self.top_right.index

    @property
    def bottom_right_index(self):
        return self.bottom_right.index

    @property
    def bottom_left_index(self):
        return self.bottom_left.index

    def __eq__(self, other):
        return self.top_left == other.top_left and self.top_right == other.top_right and self.bottom_right == other.bottom_right and self.bottom_left == other.bottom_left

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        return "Element({}, {}, {}, {})".format(repr(self.top_left),
                                                repr(self.top_right),
                                                repr(self.bottom_right),
                                                repr(self.bottom_left))

    def __str__(self):
        return "Element(tl={},tr={},br={},bl={})".format(self.top_left_index, self.top_right_index, self.bottom_right_index, self.bottom_left_index)


class Mesh(object):
    """docstring for Mesh"""

    def __init__(self, elements, nodes, material_to_element, fixed_nodes):
        self.elements = elements
        self.nodes = nodes
        self.material_to_elements = material_to_element
        self.fixed_nodes = fixed_nodes

    def nodes_count(self):
        return len(self.nodes)

    def material_for_element(self, element):
        return self.material_to_elements.get(element)

    def get_fixed_nodes_indicies(self):
        return [node.index for node in self.fixed_nodes]

    @staticmethod
    def generate(width, layers, elements_width, elements_height_per_layer, boundary_conditions):
        d_x = width / elements_width

        elements = set()
        nodes = set()
        material_to_elements = {}

        for layer in layers:
            d_y = layer.height() / elements_height_per_layer
            y = layer.top
            for i in range(elements_height_per_layer):
                x = 0
                for j in range(elements_width):
                    top_left_index = (layer.index_from_top * elements_height_per_layer * (elements_width + 1)) + i * (elements_width + 1) + j
                    top_right_index = top_left_index + 1
                    bottom_left_index = top_left_index + elements_width + 1
                    botton_right_index = bottom_left_index + 1

                    top_left = MeshNode(x, y, top_left_index)
                    top_right = MeshNode(x + d_x, y, top_right_index)
                    bottom_right = MeshNode(x + d_x, y - d_y, botton_right_index)
                    bottom_left = MeshNode(x, y - d_y, bottom_left_index)
                    nodes.add(top_left)
                    nodes.add(top_right)
                    nodes.add(bottom_right)
                    nodes.add(bottom_left)

                    element = MeshElement(top_left, top_right, bottom_right, bottom_left)
                    elements.add(element)

                    material_to_elements[element] = layer.material
                    x += d_x
                y -= d_y

        fixed_nodes_indicies = []
        if (boundary_conditions == model.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS):
            fixed_nodes_indicies = [len(nodes) - 1, len(nodes) - elements_width - 1]

        fixed_nodes = [node for node in nodes if node.index in fixed_nodes_indicies]
        return Mesh(elements, nodes, material_to_elements, fixed_nodes)
