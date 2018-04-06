class Node2D(object):
    """docstring for MeshNode"""

    def __init__(self, x1, x2, index):
        self.x1 = x1
        self.x2 = x2
        self.index = index

    def __eq__(self, other):
        return self.index == other.index

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        return "Node('{} - ({:f};{:f})')".format(self.index, self.x1, self.x2)


class FiniteElement2D(object):
    """docstring for MeshElement"""

    def __init__(self, top_left, top_right, bottom_right, bottom_left, material):
        self.top_left = top_left
        self.top_right = top_right
        self.bottom_right = bottom_right
        self.bottom_left = bottom_left
        self.material = material

    def height(self):
        return self.top_right.x2 - self.bottom_left.x2

    def width(self):
        return self.top_right.x1 - self.bottom_left.x1

    def contains(self, x1, x2):
#        if (self.top_right_index == 10):
#            print(self.top_right.x1)
        return self.bottom_left.x1 <= x1 and x1 <= self.top_right.x1 and self.bottom_left.x2 <= x2 and x2 <= self.top_right.x2

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

    def to_model_coordinates(self, ksi, teta):
        x1 = self.width() * ksi / 2 + (self.top_left.x1 + self.top_right.x1) / 2
        x2 = self.height() * teta / 2 + (self.top_left.x2 + self.bottom_left.x2) / 2
        return x1, x2

    def to_element_coordinates(self, x1, x2):
        ksi = 2 * x1 / self.width() - (self.top_left.x1 + self.top_right.x1) / self.width()
        teta = 2 * x2 / self.height() - (self.top_left.x2 + self.bottom_left.x2) / self.height()
        return ksi, teta

    def jacobian_element_coordinates(self):
        return 0.25 * self.width() * self.height()
