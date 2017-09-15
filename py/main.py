import fem.model 
import fem.mesh
import fem.solver


# print(dir(fem.model))
# print(fem.model.Model)

l = 1
K = 0.8
h = 0.1

gA=0.03
gV=20

N = 10
M = 2

geometry = fem.model.Geometry(l,K, gA, gV)
layer = fem.model.Layer(-h/2, h/2, fem.model.Material.steel())
layers = [layer]
print(layers)
model = fem.model.Model(geometry, layers, fem.model.Model.SHARNIR)

mesh = fem.mesh.Mesh.generate(model.geometry.width, layers, N, M)
print(mesh.elements)
print(len(mesh.elements))


# [vec lam] = solve(geom, layerModel, N, M, staticIndecies);

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
# x=0:l/N:l;
# h = plot(x, y);
# axis([0 l -1 1]);
# for t = 0.001:0.0001:5
#   y = (cos(w*t)+sin(w*t)).*midPaneResult;
#   set(h, 'YData', y);
#   pause(0.1);
# end