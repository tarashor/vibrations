function [ oArgs ] = MatrixG(l, K, gA, gV, alpha1, alpha2)

N=zeros(2,2);
N(1,1) = g11(l, K, gA, gV, alpha1, alpha2);
N(2,2) = 1;

oArgs = N;

end
