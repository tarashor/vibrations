function [ oArgs ] = localStiffnessMatrix(psi, teta, E, v, K, alpha1start, alpha1end, alpha2start, alpha2end)

alpha2 = (alpha2end-alpha2start)*teta / 2 + (alpha2end+alpha2start) / 2;
A=MatrixA(E,v,K,alpha2);
B=MatrixB(K,alpha2);
I=MatrixI(alpha1start, alpha1end, alpha2start, alpha2end);
N=MatrixN(psi, teta);
J=0.25*(alpha1end-alpha1start)*(alpha2end-alpha2start);


lSM = N'*I'*B'*A*B*I*N*J;
oArgs = lSM;

end
