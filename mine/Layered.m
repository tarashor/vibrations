clc;
clear;

N = 20;
l = 1;
K = 0;

h1 = 0.1;
E1 = 100000;
v1=0.3;
rho1=9000;

h2 = 0.3;
E2 = 100;
v2=0.5;
rho2=7000;

layers(1) = GetModel(0,h1,l,K,rho1,E1,v1,N);
layers(2) = GetModel(h1,h1+h2,l,K,rho2,E2,v2,N);
layers(3) = GetModel(h1+h2,h1+h2+h1,l,K,rho1,E1,v1,N);

[meshBegin, meshEnd] = GenerateMesh(N, l);
x = [meshBegin meshEnd(N)];

for k=1:length(layers)
  s = StiffnessMatrix(layers(i), meshBegin, meshEnd);
  m = rho*MassMatrix(layers(i), meshBegin, meshEnd);
end

staticIndecies = getBoundaryConditionIndicies(N);
s = applyStaticBoundaryConditionsToMatrix(s, staticIndecies);
m = applyStaticBoundaryConditionsToMatrix(m, staticIndecies);

 [vec, lam] = eig (s, m);
  [vec, lam] = sortResults(vec, lam);
  
  eigvec = extendresultWithStaticBoundaryConditions(vec, staticIndecies);
  eigval = lam;