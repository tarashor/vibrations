clc;
clear;

N = 20;
l = 1;
K = 0;

geom=[l,K];

h1 = 0.02;
[E1 v1 rho1] = GetSteel();

h2 = 0.01;
[E2 v2 rho2] = GetRubber();

layersCount=3;
layers=cell(1,layersCount);
layers{1} = GetLayerModel(0, h1, rho1, E1, v1);
layers{2} = GetLayerModel(h1, h1+h2, rho2, E2, v2);
layers{3} = GetLayerModel(h1+h2, h1+h2+h1, rho1, E1, v1);

staticIndecies = getBoundaryConditionIndiciesForLayeredMatrix(N, layersCount);
[vec lam x] = solveLayered(geom, layers, N, staticIndecies);

ind = 1;
lam(ind)
resVector = vec(:, ind);

%layerToShow = fix(K./2)+1;
%midPaneResult=zeros(N+1,1);
%for l=1:N+1
%  i_new = (2*K+1)*(2*l-1)+(layerToShow-1)*2;
%  temp=(resVector(i_new+1)+resVector(i_new+3))./2+resVector(i_new+2);
%  midPaneResult(l) = temp;
%end

layerToShow = 1;
midPaneResult=zeros(N+1,1);
for p=1:N+1
  i_new = (2*layersCount+1)*(2*p-1)+(layerToShow-1)*2;
  temp=resVector(i_new+1);
  midPaneResult(p) = temp;
end

%figure(ind);
%plot(x, midPaneResult);

w=sqrt(lam(ind));
y = cos(w*0).*midPaneResult;
h = plot(x, y);
axis([0 l -1 1]);
for t = 0.001:0.0001:5
  y = cos(w*t).*midPaneResult;
  set(h, 'YData', y);
  pause(0.1);
end


