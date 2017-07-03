function [ oArgs ] = MatrixA(E,v,K,alpha2)
q=1+K*alpha2;
q4=q^4;
q2=q^2;

A=zeros(3,3);


A(1,1)=E*(1-v)/(q4*(1+v)*(1-2*v));
A(2,2)=E*(1-v)/((1+v)*(1-2*v));
A(1,2)=E*v/(q2*(1+v)*(1-2*v));
A(2,1)=E*v/(q2*(1+v)*(1-2*v));
A(3,3)=E/(q2*(1+v)*2);

oArgs = A;

end