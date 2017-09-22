import fem.model
import fem.mesh
import fem.solver
import matplotlib.pyplot as plt


# print(dir(fem.model))
# print(fem.model.Model)

width = 1
curvature = 0.8
thickness = 0.1

corrugation_amplitude = 0.03
corrugation_frequency = 20
# corrugation_amplitude = 0
# corrugation_frequency = 0

layers_count = 1

N = 40
M = 11

geometry = fem.model.Geometry(width, curvature, corrugation_amplitude, corrugation_frequency)

layer_top = thickness / 2
layer_thickness = thickness / layers_count
layers = set()
for i in range(layers_count):
    layer = fem.model.Layer(layer_top - layer_thickness, layer_top, fem.model.Material.steel(), i)
    layers.add(layer)
    layer_top -= layer_thickness

# layers = tuple([layer])

# print(layers)

model = fem.model.Model(geometry, layers, fem.model.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)

mesh = fem.mesh.Mesh.generate(model.geometry.width, layers, N, M, model.boundary_conditions)

# list_nodes = sorted(mesh.nodes, key=lambda n: n.index)
# for n in list_nodes:
#    print(n)
#
# list__fixed_nodes = sorted(mesh.get_fixed_nodes_indicies())
# for i in list__fixed_nodes:
#    print(i)


result = fem.solver.solve(model, mesh)
l, v1, v2 = result.get_result(0)

print(l)

x = []
y = []

list_nodes = sorted(mesh.nodes, key=lambda n: n.index)
print("Nodes = {}".format(len(mesh.nodes)))

fn = mesh.fixed_nodes
for n in fn:
    print(n.x + v1[n.index])
    print(n.y + v2[n.index])

for n in list_nodes:
    x.append(n.x + v1[n.index])
    y.append(n.y + v2[n.index])

plt.plot(x, y, 'ro')
# plt.axis([-2, 2, -2, 2])
plt.show()

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
# x=0:width/N:width;
# height = plot(x, y);
# axis([0 width -1 1]);
# for t = 0.001:0.0001:5
#   y = (cos(w*t)+sin(w*t)).*midPaneResult;
#   set(height, 'YData', y);
#   pause(0.1);
# end
