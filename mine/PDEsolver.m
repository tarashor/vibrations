clc;
clear;

N = 20;
l = 1;
h = 0.1;
K = 0;
E = 100000;
v=0.3;
rho=9000;

iArgs = GetModel(0,h,l,K,rho,E,v,N);
staticIndecies = getBoundaryConditionIndicies(N);
[vec lam xVector] = solve(iArgs, staticIndecies);


%count = length(lam);
%countEnd = count-2;
%for j=count:-1:countEnd
%  ind=j;
  ind = 1;
  lam(ind)
  resVector = vec(:, ind);
  for i=1:N+1
      z(i)=x(resVector, xVector, [xVector(i),0,-h/2], h);
  end
  figure(ind);
  plot(xVector, z);
%end