clc;
clear;

N = 20;
l = 1;
h = 0.01;
E = 100000;
v=0.3;
rho=8000;

staticIndecies = getBoundaryConditionIndicies(N);
KMax = 1;
Ks = 0:0.05:KMax;
lMin = zeros(length(Ks),1);
for i=1:length(Ks)
  iArgs = GetModel(h,l,Ks(i),rho,E,v,N);
  [vec lam xVector] = solve(iArgs, staticIndecies);
  lMin(i) = lam(1);
end

plot(Ks, lMin);