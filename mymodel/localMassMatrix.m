function [ oArgs ] = localMassMatrix(psi, teta, l, K, gA, gV, alpha1start, alpha1end, alpha2start, alpha2end)

  N2=MatrixN2(psi, teta);
  
  alpha1 = (alpha1end-alpha1start)*teta / 2 + (alpha1end+alpha1start) / 2;
  alpha2 = (alpha2end-alpha2start)*teta / 2 + (alpha2end+alpha2start) / 2;
  G=MatrixG(l, K, gA, gV, alpha1, alpha2);
  
  J=0.25*(alpha1end-alpha1start)*(alpha2end-alpha2start);

  lMM = N2'*G*N2*J;
  
  oArgs = lMM;

end
