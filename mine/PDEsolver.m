clc;
clear;

N = 20;
l = 1;
K = 0;

h = 0.04;
E = 210000000000;
v=0.3;
rho=8000;

iArgs = GetModel(0,h,l,K,rho,E,v,N);
staticIndecies = getBoundaryConditionIndicies(N);
[vec lam xVector] = solve(iArgs, staticIndecies);


%count = length(lam);
%countEnd = count-2;
%for j=count:-1:countEnd
%  ind=j;
  ind = 4;
  lam(ind)
  resVector = vec(:, ind);
  for i=1:N+1
      z(i)=x(resVector, xVector, [xVector(i),0,-h/2], h);
  end
  %figure(ind);
  %plot(xVector, z);
  
  w=sqrt(lam(ind));
  y = cos(w*0).*z;
  h = plot (xVector, y);
  axis([0 l -1 1]);
  for t = 0.001:0.0001:5
    y = cos(w*t).*z;
    set (h, 'YData', y);
    pause(0.1);
  endfor
  
%end