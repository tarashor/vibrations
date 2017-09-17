function [ oArgs ] = localStiffnessMatrix(psi, teta, l, K, gA, gV, E, v, alpha1start, alpha1end, alpha2start, alpha2end)

alpha1 = (alpha1end-alpha1start)*psi / 2 + (alpha1end+alpha1start) / 2;
alpha2 = (alpha2end-alpha2start)*teta / 2 + (alpha2end+alpha2start) / 2;
A=MatrixA(E,v, l, K, gA, gV, alpha1, alpha2);
%T=MatrixT(l, K, gA, gV, alpha1, alpha2);
B=MatrixB(l, K, gA, gV, alpha1, alpha2);

I=MatrixI(alpha1start, alpha1end, alpha2start, alpha2end);
N=MatrixN(psi, teta);
J=0.25*(alpha1end-alpha1start)*(alpha2end-alpha2start);

lSM = N'*I'*B'*A*B*I*N*J;
oArgs = lSM;

end
