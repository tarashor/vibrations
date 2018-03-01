clear;
clc;
g_a = 0.03;
g_f = 20;
K=0.8;
L = 2;
h=0.05;

num_points = 200;

R = 1/K;
delta = L/num_points;
a1 = 0:delta:L;

a2 = [-h/2, 0, h/2];

ar = pi/2+(L-2*a1)/(2*R);

r(1,:)=(R + g_a.*cos(g_f.*ar)).*cos(ar);
r(2,:)=(R + g_a.*cos(g_f.*ar)).*sin(ar);

a=g_a*g_f*K.*sin(g_f.*ar);
b=(1+g_a*K.*cos(g_f.*ar));

n(1,:)=-a.*sin(ar)+b.*cos(ar);
n(2,:)=a.*cos(ar)+b.*sin(ar);

length_n=sqrt(a.*a+b.*b); 

n=n./length_n;

l = sqrt(n(1,:).*n(1,:)+n(2,:).*n(2,:));

figure;
hold on;  
for i=1:length(a2)
  R=r+a2(i)*n;
  plot(R(1,:),R(2,:), 'g');
end

alpha2 = a2(2);

x=r(1,:)+alpha2*n(1,:);
y=r(2,:)+alpha2*n(2,:);

l2 = a.*a+b.*b;
dl = sqrt(l2.*l2.*l2);
%
%dn(1,:)=l2.*(2*g_a*g_f*K*K.*sin(g_f*ar).*cos(ar)+(g_a*K*K*(1+g_f*g_f).*cos(g_f*ar)+K).*sin(ar));
dn(2,:)=l2.*(2*g_a*g_f*K*K.*sin(g_f*ar).*sin(ar)-(g_a*K*K*(1+g_f*g_f).*cos(g_f*ar)+K).*cos(ar));
%dn(1,:)=dn(1,:)-(1+g_a*(1-g_f*g_f)*K.*cos(g_f*ar)).*(g_a*g_f*K*K.*sin(g_f*ar)).*((1+g_a*K.*cos(g_f.*ar)).*cos(ar)-g_a*g_f*K*sin(g_f.*ar).*sin(ar));
dn(2,:)=dn(2,:)-(1+g_a*(1-g_f*g_f)*K.*cos(g_f*ar)).*(g_a*g_f*K*K.*sin(g_f*ar)).*((1+g_a*K.*cos(g_f.*ar)).*sin(ar)+g_a*g_f*K*sin(g_f.*ar).*cos(ar));


dn(1,:)=g_a*g_f*K*K.*sin(g_f*ar).*cos(ar).*(1+2*g_f*g_f*g_a*g_a*K*K+(2+g_f*g_f)*g_a*K.*cos(g_f*ar)+g_a*g_a*K*K.*cos(g_f*ar).*cos(g_f*ar).*(1-g_f*g_f));
dn(1,:)=dn(1,:)+sin(ar).*(K+g_a*K*K.*cos(g_f*ar).*(3+g_f*g_f+2*g_f*g_f*g_f*g_f*g_a*g_a*K*K)+g_a*g_a*K*K*K.*cos(g_f*ar).*cos(g_f*ar)*(3+2*g_f*g_f)+g_a*g_a*g_a*K*K*K*K.*cos(g_f*ar).*cos(g_f*ar).*cos(g_f*ar)*(1+g_f*g_f-2*g_f*g_f*g_f*g_f));

%dn(2,:)=g_a*g_f*K*K.*sin(g_f*ar).*cos(ar).*(1+2*g_f*g_f*g_a*g_a*K*K+(2+g_f*g_f)*g_a*K.*cos(g_f*ar)+g_a*g_a*K*K.*cos(g_f*ar).*cos(g_f*ar).*(1-g_f*g_f));
%dn(2,:)=dn(1,:)-sin(ar).*(K+g_a*K*K.*cos(g_f*ar).*(3+g_f*g_f+2*g_f*g_f*g_f*g_f*g_a*g_a*K*K)+g_a*g_a*K*K*K.*cos(g_f*ar).*cos(g_f*ar)*(3+2*g_f*g_f)+g_a*g_a*g_a*K*K*K*K.*cos(g_f*ar).*cos(g_f*ar).*cos(g_f*ar)*(1+g_f*g_f-2*g_f*g_f*g_f*g_f));


dn=dn./dl;

dr(1,:)=a.*cos(ar)+b.*sin(ar);
dr(2,:)=a.*sin(ar)-b.*cos(ar);

R1=dr+alpha2*dn;
R2=n;

ip=3;

%quiver (x(1:ip:end),y(1:ip:end), v(1, 1:ip:end), v(2, 1:ip:end), 'k')

%quiver (x(1:ip:end),y(1:ip:end), dr(1, 1:ip:end), dr(2, 1:ip:end), 'm')
%quiver (x(1:ip:end),y(1:ip:end), dn(1, 1:ip:end), dn(2, 1:ip:end), 'c')

%R1'*R2

l = R1(1,:).*R2(1,:)+R1(2,:).*R2(2,:);

%quiver (x(1:ip:end),y(1:ip:end), R1(1, 1:ip:end), R1(2, 1:ip:end), 'r')
quiver (x(1:ip:end),y(1:ip:end), R2(1, 1:ip:end), R2(2, 1:ip:end), 'b')

daspect([1 1 1])
