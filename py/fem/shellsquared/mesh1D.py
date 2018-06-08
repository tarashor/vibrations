from . import finiteelements1D as fe

class Mesh1D(object):
    """docstring for Mesh"""

    def __init__(self, elements, nodes, fixed_nodes):
        self.elements = elements
        self.nodes = nodes
        self.fixed_nodes = fixed_nodes

    def nodes_count(self):
        return len(self.nodes)

    def get_element(self, x1):
        element = None
        for el in self.elements:
            if (el.contains(x1)):
                element = el
                break

        return element

    def get_fixed_nodes_indicies(self):
        return [node.index for node in self.fixed_nodes]

    @staticmethod
    def generate(width, layers, elements_width, boundary_conditions):
        d_x = width / elements_width

        elements = set()
        nodes = set()

        for layer in layers:
            x = 0
            for j in range(elements_width):
                start_index = (layer.index_from_top * (elements_width + 1)) + j
                end_index = start_index + 1

                start = fe.Node1D(x, start_index)
                end = fe.Node1D(x + d_x, end_index)
                nodes.add(start)
                nodes.add(end)

                element = fe.FiniteElement1D(start, end, layer.material)
                elements.add(element)
                x += d_x
        
        fixed_nodes_indicies = []
        if (boundary_conditions == 1):
            fixed_nodes_indicies = [0, len(nodes) - 1]

        fixed_nodes = [node for node in nodes if node.index in fixed_nodes_indicies]
        return Mesh1D(elements, nodes, fixed_nodes)
