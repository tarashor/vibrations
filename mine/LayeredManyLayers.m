clc;
clear;

R = 2;
fi0 = pi/4;
layersCount=11;
lDivHCoef = 0.1;

l = 2*fi0*R
hTotal = l*lDivHCoef;
K = 1/R;
hLayer = hTotal/layersCount;
N = 20;
[ESteel vSteel rhoSteel] = GetSteel();
[ERubber vRubber rhoRubber] = GetRubber();

geom=[l,K];

layers=cell(1,layersCount);
hCurrent = 0;
for i=1:layersCount
  E = ESteel;
  v = vSteel;
  rho = rhoSteel;
  if (rem(i,2) == 0)
    E = ERubber;
    v = vRubber;
    rho = rhoRubber;
  endif
  layers{i} = GetLayerModel(hCurrent, hCurrent + hLayer, rho, E, v);  
  hCurrent += hLayer
end

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
midPaneResult
y = cos(w*0).*midPaneResult;
h = plot(x, y);
axis([0 l -1 1]);
for t = 0.001:0.0001:5
  y = cos(w*t).*midPaneResult;
  set(h, 'YData', y);
  pause(0.01);
end