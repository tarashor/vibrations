function [ oArgs ] = MatrixB(K,alpha2)
q=1+K*alpha2;

B=zeros(3,6);
B(1,2)=B(2,6)=B(3,3)=B(3,5)=1;
B(1,4)=K*q;
B(3,1)=-2*K/q;
oArgs = B;

end