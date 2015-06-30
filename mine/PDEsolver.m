clc;
clear;
l = 1;
h = 0.01;
K = 0;
E = 100000;
v=0.3;
N=20;
rho=8000;
iArgs = [h,l,K,E,v,N];
[meshBegin, meshEnd] = GenerateMesh(N, l);
xVector=[meshBegin meshEnd(N)];

s=StiffnessMatrix(iArgs, meshBegin, meshEnd);
m = rho*MassMatrix(iArgs, meshBegin, meshEnd);

s = applyStaticBoundaryConditionsToMatrix(s, N);
m = applyStaticBoundaryConditionsToMatrix(m, N);

[vec, lam] = eig (s, m);
[vec, lam] = sortResults(vec, lam);
ind=3;

%lam

%s*vec(:, ind)
%lam(ind,ind)*m*vec(:, ind)

xVector = [meshBegin meshEnd(N)];
resVector = extendresultWithStaticBoundaryConditions(vec(:, ind), N);

resVector=resVector/norm(resVector);
%q=[00.000   01.465   -00.024   00.000   00.005   -00.030   00.000   00.000   00.000   -05.443   -05.355   -00.008   00.000   -01.465   00.024   00.000   00.005   -00.030];
%q=q/norm(q)

for i=1:N+1
    z(i)=x(resVector, xVector, [xVector(i),0,-h/2], h);
end
z

plot(xVector, z);