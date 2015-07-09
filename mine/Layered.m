clc;
clear;

N = 20;
l = 1;
K = 1;

h1 = 0.02;
E1 = 210000000000;
v1=0.3;
rho1=8000;

h2 = 0.01;
E2 = 10000000;
v2=0.49;
rho2=1150;

K=3;

layers=cell(1,K);

layers{1} = GetModel(0,h1,l,K,rho1,E1,v1,N);
layers{2} = GetModel(h1,h1+h2,l,K,rho2,E2,v2,N);
layers{3} = GetModel(h1+h2,h1+h2+h1,l,K,rho1,E1,v1,N);

dim = 2*(N+1)*(2*K+1);

[meshBegin, meshEnd] = GenerateMesh(N, l);
x = [meshBegin meshEnd(N)];

SMatrix = zeros(dim,dim);
MMatrix = zeros(dim,dim);

for k=1:K
  s = StiffnessMatrix(layers{k}, meshBegin, meshEnd);
  m = MassMatrix(layers{k}, meshBegin, meshEnd);
  
  SMatrix = SumLayerMatrix(SMatrix, s, k, K);
  MMatrix = SumLayerMatrix(MMatrix, m, k, K);
end

staticIndecies = getBoundaryConditionIndiciesForLayeredMatrix(N, K);
SMatrix = applyStaticBoundaryConditionsToMatrix(SMatrix, staticIndecies);
MMatrix = applyStaticBoundaryConditionsToMatrix(MMatrix, staticIndecies);

[vec, lam] = eig (SMatrix,MMatrix);
[vec, lam] = sortResults(vec, lam);

eigvec = extendresultWithStaticBoundaryConditions(vec, staticIndecies);

ind = 1;
lam(ind)
resVector = eigvec(:, ind);

%layerToShow = fix(K./2)+1;
%midPaneResult=zeros(N+1,1);
%for l=1:N+1
%  i_new = (2*K+1)*(2*l-1)+(layerToShow-1)*2;
%  temp=(resVector(i_new+1)+resVector(i_new+3))./2+resVector(i_new+2);
%  midPaneResult(l) = temp;
%end

layerToShow = 1;
midPaneResult=zeros(N+1,1);
for l=1:N+1
  i_new = (2*K+1)*(2*l-1)+(layerToShow-1)*2;
  temp=resVector(i_new+1);
  midPaneResult(l) = temp;
end

%figure(ind);
%plot(x, midPaneResult);

w=sqrt(lam(ind));
  y = cos(w*0).*midPaneResult;
  h = plot (x, y);
  axis([0 l -1 1]);
  for t = 0.001:0.0001:5
    y = cos(w*t).*midPaneResult;
    set (h, 'YData', y);
    pause(0.1);
  endfor


