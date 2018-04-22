function [ oArgs ] = MatrixA(E,v,l, K, gA, gV, alpha1, alpha2)

A=zeros(3,3);

g=g11(l, K, gA, gV, alpha1, alpha2);

A(1,1)=E*(1-v)*g*g/((1+v)*(1-2*v));
A(2,2)=E*(1-v)/((1+v)*(1-2*v));
A(1,2)=E*v*g/((1+v)*(1-2*v));
A(2,1)=E*v*g/((1+v)*(1-2*v));
A(3,3)=E*g/((1+v)*2);

oArgs = A;

end