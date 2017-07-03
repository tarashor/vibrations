function [ oArgs ] = MatrixG(K, alpha2)

q=1+K*alpha2;
q2=q^2;

N=zeros(2,2);
N(1,1) = 1/q2;
N(2,2) = 1;

oArgs = N;

end
