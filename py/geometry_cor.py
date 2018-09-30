import fem.geometry as g


import plot

width = 2
curvature = 0.8
thickness = 0.05

corrugation_amplitude = 0
corrugation_frequency = 20

geometry1 = g.General(width, curvature, corrugation_amplitude, corrugation_frequency)

geometry2 = g.General(width, curvature, 0, 0)

#plot.plot_midpane_in_cartesian(0, width, -thickness / 2, thickness / 2, geometry1.to_cartesian_coordinates, geometry2.to_cartesian_coordinates)

plot.plot_init_geometry(geometry1, 0, width, -thickness / 2, thickness / 2, 0)