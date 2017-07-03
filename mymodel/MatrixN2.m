function [ oArgs ] = MatrixN2(psi, teta)

N=zeros(2,8);
N(1,1) = N(2,5) = 0.25*(1-psi)*(1+teta);
N(1,2) = N(2,6) = 0.25*(1+psi)*(1+teta);
N(1,3) = N(2,7) = 0.25*(1+psi)*(1-teta);
N(1,4) = N(2,8) = 0.25*(1-psi)*(1-teta);

oArgs = N;

end
