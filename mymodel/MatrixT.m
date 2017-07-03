function [ oArgs ] = MatrixT(l, K, gA, gV, alpha1, alpha2)

T=zeros(6,6);
g=g11(l, K, gA, gV, alpha1, alpha2);
T(1,1) = g*g;
T(2,2) = 1;
T(3,3) = 1;
T(4,4) = g;
T(5,5) = 1;
T(6,6) = g;

oArgs = T;

end
