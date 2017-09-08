
class Geometry:
	def __init__(self, length, curvature, wave_amp, wave_frequency):
		self.length = length
		self.curvature = curvature
		self.wave_amp = wave_amp
		self.wave_frequency = wave_frequency

class Material:
	def __init__(self, E, v, rho):
		self.E = E
		self.v = v
		self.rho = rho

	@staticmethod
	def steel():
		return Material(210000000000, 0.3, 8000)

class Layer:
	def __init__(self, bottomX, topX, material):
		self.bottomX = bottomX
		self.topX = topX
		self.material = material	


l = 2
K = 0.8
h = 0.05

gA=0.03
gV=20

N = 40
M = 6

geometry = Geometry(l,K, gA, gV)
material = Material.steel()
layer = Layer(-h/2, h/2, material)

# staticIndecies = getBoundaryConditionIndicies(N, M);
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