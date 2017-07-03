function [ oArgs ] = MatrixB(l, K, gA, gV, alpha1, alpha2)
q=1+K*alpha2;
a=(pi+K*l)/2-K*alpha1;
w = q+gA*K*cos(gV*a);
z = gA*gV*K*sin(gV*a);

G111=2*z*K/w;
G211=-gA*gV*gV*K*K*cos(gV*a)-w*K-2*z*z*K/w;
G112=K/q;

B=zeros(3,6);
B(1,1)=-G111;
B(1,4)=-G211;
B(1,2)=B(2,6)=B(3,3)=B(3,5)=1;
B(3,1)=-2*G112;
oArgs = B;

end