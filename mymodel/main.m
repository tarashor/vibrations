clc;
clear;

l = 2;
K = 0.8;
h= 0.05;

N = 40;
M=6;

geom=[l,K];
[E v rho] = GetSteel();

layerModel = GetLayerModel(-h/2, h/2, rho, E, v);

staticIndecies = getBoundaryConditionIndicies(N, M);
[vec lam] = solve(geom, layerModel, N, M, staticIndecies);

ind = 1;

printf ("Minimum frequancy = %f\n", sqrt(lam(ind)/rho));
resVector = vec(:, ind);

printf ("Norm of resVector = %f\n", sqrt(resVector'*resVector));

%layerToShow = fix(K./2)+1;
%midPaneResult=zeros(N+1,1);
%for l=1:N+1
%  i_new = (2*K+1)*(2*l-1)+(layerToShow-1)*2;
%  temp=(resVector(i_new+1)+resVector(i_new+3))./2+resVector(i_new+2);
%  midPaneResult(l) = temp;
%end

midPaneResult=zeros(N+1,1);
for p=1:N+1
  i_new = (M+1)*(N+1)+(M/2)*(N + 1) + p;
  temp=resVector(i_new);
  midPaneResult(p) = temp;
end

%figure(ind);
%plot(x, midPaneResult);

w=sqrt(lam(ind)/rho);
y = (cos(w*0)+sin(w*0)).*midPaneResult;
x=0:l/N:l;
h = plot(x, y);
axis([0 l -1 1]);
for t = 0.001:0.0001:5
  y = (cos(w*t)+sin(w*t)).*midPaneResult;
  set(h, 'YData', y);
  pause(0.1);
end
