class Node1D(object):
    """docstring for MeshNode"""

    def __init__(self, x1, index):
        self.x1 = x1
        self.index = index

    def __eq__(self, other):
        return self.index == other.index

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        return "Node1D('{} - ({:f})')".format(self.index, self.x1)


class FiniteElement1D(object):
    
    """docstring for MeshElement"""

    def __init__(self, start, end, material):
        self.start = start
        self.end = end
        self.material = material

    def width(self):
        return self.end.x1 - self.start.x1

    def contains(self, x1):
        eps = 0.00000001
#        if (self.top_right_index == 10):
#            print(self.top_right.x1)
        return self.start.x1 - eps <= x1 and x1 <= self.end.x1 + eps

    @property
    def start_index(self):
        return self.start.index

    @property
    def end_index(self):
        return self.end.index

    def __eq__(self, other):
        return self.start == other.start and self.end == other.end

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        return "Element1D({}, {})".format(repr(self.start),
                                                repr(self.end))

    def __str__(self):
        return "Element1D(start={},end={})".format(self.start_index, self.end_index)

    def to_model_coordinates(self, ksi):
        x1 = self.width() * ksi / 2 + (self.start.x1 + self.end.x1) / 2
        return x1

    def to_element_coordinates(self, x1):
        ksi = 2 * x1 / self.width() - (self.start.x1 + self.end.x1) / self.width()
        return ksi

    def jacobian_element_coordinates(self):
        return 0.5 * self.width()
