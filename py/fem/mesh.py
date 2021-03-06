from . import model
from . import finiteelements as fe


class Mesh(object):
    """docstring for Mesh"""

    def __init__(self, elements, nodes, fixed_nodes):
        self.elements = elements
        self.nodes = nodes
        self.fixed_nodes = fixed_nodes

    def nodes_count(self):
        return len(self.nodes)

    def get_element(self, x1, x2):
        element = None
        for el in self.elements:
            if (el.contains(x1, x2)):
                element = el
                break

        return element

    def get_fixed_nodes_indicies(self):
        return [node.index for node in self.fixed_nodes]

    @staticmethod
    def generate(width, layers, elements_width, elements_height_per_layer, boundary_conditions):
        d_x = width / elements_width

        elements = set()
        nodes = set()

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

                    top_left = fe.Node2D(x, y, top_left_index)
                    top_right = fe.Node2D(x + d_x, y, top_right_index)
                    bottom_right = fe.Node2D(x + d_x, y - d_y, botton_right_index)
                    bottom_left = fe.Node2D(x, y - d_y, bottom_left_index)
                    nodes.add(top_left)
                    nodes.add(top_right)
                    nodes.add(bottom_right)
                    nodes.add(bottom_left)

                    element = fe.FiniteElement2D(top_left, top_right, bottom_right, bottom_left, layer.material)
                    elements.add(element)

                    x += d_x
                y -= d_y

        fixed_nodes_indicies = []
        if (boundary_conditions == model.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS):
            fixed_nodes_indicies = [len(nodes) - 1, len(nodes) - elements_width - 1]

        fixed_nodes = [node for node in nodes if node.index in fixed_nodes_indicies]
        return Mesh(elements, nodes, fixed_nodes)
