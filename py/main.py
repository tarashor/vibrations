import fem.model
import fem.mesh
import fem.solver


# print(dir(fem.model))
# print(fem.model.Model)

width = 1
curvature = 0.8
thickness = 0.1

corrugation_amplitude = 0.03
corrugation_frequency = 20

layers_count = 1

N = 10
M = 2

geometry = fem.model.Geometry(width, curvature, corrugation_amplitude, corrugation_frequency)

layer_top = thickness / 2
layer_thickness = thickness / layers_count
layers = set()
for i in range(layers_count):
    layer = fem.model.Layer(layer_top - layer_thickness, layer_top, fem.model.Material.steel())
    layers.add(layer)
    layer_top -= layer_thickness

# layers = tuple([layer])

print(layers)
model = fem.model.Model(geometry, layers, fem.model.Model.FIXED_BOTTOM_LEFT_RIGHT_POINTS)

mesh = fem.mesh.Mesh.generate(model.geometry.width, layers, N, M)
for element in mesh.nodes:
    print(str(element))
print(len(mesh.elements))


result = fem.solver.solve(model, mesh)
print(result.get_result(0))

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
