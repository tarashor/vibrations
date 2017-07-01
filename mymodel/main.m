clc;
clear;

l = 1;
K = 0;
h= 0.04;

N = 20;
M=6;

geom=[l,K];
[E v rho] = GetSteel();

layerModel = GetLayerModel(-h/2, h/2, rho, E, v);

staticIndecies = getBoundaryConditionIndiciesForLayeredMatrix(N, M);
[vec lam x] = solve(geom, layerModel, N, M, staticIndecies);

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
