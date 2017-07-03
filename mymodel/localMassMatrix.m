function [ oArgs ] = localMassMatrix(psi, teta, K, alpha1start, alpha1end, alpha2start, alpha2end)

alpha2 = (alpha2end-alpha2start)*teta / 2 + (alpha2end+alpha2start) / 2;
N2=MatrixN2(psi, teta);
G=MatrixG(K, alpha2);
J=0.25*(alpha1end-alpha1start)*(alpha2end-alpha2start);


lMM = N2'*G*N2*J;
oArgs = lMM;

end
