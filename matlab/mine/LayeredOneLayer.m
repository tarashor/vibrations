clc;
clear;


l = 2;
K = 0.8;
h= 0.05;

N = 40;

geom=[l,K];
h1 = h;
[E1 v1 rho1] = GetSteel();

layersCount=1;
layers=cell(1,layersCount);
layers{1} = GetLayerModel(0, h1, rho1, E1, v1);

staticIndecies = getBoundaryConditionIndiciesForLayeredMatrix(N, layersCount);
[vec lam x] = solveLayered(geom, layers, N, staticIndecies);

ind = 1;
printf ("Minimum frequancy = %f\n", sqrt(lam(ind)/rho1));
resVector = vec(:, ind);

printf ("Norm of resVector = %f\n", sqrt(resVector'*resVector));

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

w=sqrt(lam(ind)/rho1);
y = (cos(w*0)+sin(w*0)).*midPaneResult;
h = plot(x, y);
axis([0 l -1 1]);
for t = 0.001:0.0001:5
  y = (cos(w*t)+sin(w*t)).*midPaneResult;
  set(h, 'YData', y);
  pause(0.1);
end


