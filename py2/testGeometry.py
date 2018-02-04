import fem.geometry as g
import fem.matrices as m


width = 2
curvature = 0.8
geometry = g.CylindricalPlate(width, curvature)
G = geometry.kristophel_symbols(0.09382014999999999, -0.022654496250000003, 0)
print(G)
