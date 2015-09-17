clc;
clear;

N = 20;
l = 1;
h = 0.01;
K=0;
E = 100000;
v=0.3;
rho=8000;

staticIndecies = getBoundaryConditionIndicies(N);
LMax = 1;
Ls = 0.5:0.2:2;
lMin = zeros(length(Ls),1);
for i=1:length(Ls)
  iArgs = GetModel(h,Ls(i),K,rho,E,v,N);
  [vec lam xVector] = solve(iArgs, staticIndecies);
  lMin(i) = lam(1);
end

plot(Ls, lMin);