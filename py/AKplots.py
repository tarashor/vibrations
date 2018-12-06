import fem.geometry as g
import matplotlib.pyplot as plt

width = 2
curvature = 1.25
thickness = 0.05

corrugation_amplitudes = [0, 0.015, 0.03, 0.06, 0.1, 0.2, 0.25, 0.3]
corrugation_frequencies = [2, 6, 15, 20, 50, 100]

#geometry = g.General(width, curvature, corrugation_amplitude, 24)
geometry = g.General(width, curvature, corrugation_amplitudes[7], corrugation_frequencies[2])

plot_x1_elements = 1000

x1_start = 0
x1_end = width

dx1 = (x1_end - x1_start) / plot_x1_elements

X_init = []
Y_init = []


for i in range(plot_x1_elements + 1):
    x1 = x1_start + i * dx1
    A, K = geometry.get_A_and_K(x1, 0, 0)


    X_init.append(x1)
#    Y_init.append(A)
    Y_init.append(1/(1+thickness*K))

    
plt.plot(X_init, Y_init, "r")
#plt.title("Displacement related to minimal natural frequency")
# plt.title(r"Форма панелі $L={}, h={}, K={}, g_A={}, g_v={}$".format(x1_end - x1_start, x2_end - x2_start, result.geometry.curvature, result.geometry.corrugation_amplitude, result.geometry.corrugation_frequency))
plt.axes().set_aspect('equal', 'datalim')
plt.legend(loc='best')
plt.xlabel(r"$\alpha_1$, m", fontsize=12)
plt.ylabel(r"$A$, m", fontsize=12)
plt.grid()
plt.show()
    

