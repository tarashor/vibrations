clc;
clear;

l = 2;
K = 0.8;
h= 0.05;

gA=0.03;
%gV=20;

N = 40;
M=6;

gV=[1,2,3,4,5,6,7, 8, 9, 10,12,15,20];

for i=1:length(gV)
  geom=[l,K, gA, gV(i)];
  [E v rho] = GetSteel();

  layerModel = GetLayerModel(-h/2, h/2, rho, E, v);

  staticIndecies = getBoundaryConditionIndicies(N, M);
  [vec lam] = solve(geom, layerModel, N, M, staticIndecies);

  ind = 1;
  printf ("Corrugated frequency = %d, Minimum frequency = %f\n", gV(i), sqrt(lam(ind)/rho));
end
