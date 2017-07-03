function [ oArgs ] = MatrixI(alpha1start, alpha1end, alpha2start, alpha2end)

I=eye(6);
I(2,2) = I(5,5) = 2/(alpha1end-alpha1start);
I(3,3) = I(6,6) = 2/(alpha2end-alpha2start);
oArgs = I;

end