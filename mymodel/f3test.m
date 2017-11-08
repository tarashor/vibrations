clear;
clc;
g_a = 0.03;
g_f = 20;
K=0.8;

R = 1/K;

L = 2;

p = 200;


delta = L/p;

a1 = 0:delta:L;
a2 = [-0.025, 0, 0.025];

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
for i=1:length(a2)
  z=a2(i)*n;
  %l = sqrt(z(1,:).*z(1,:)+z(2,:).*z(2,:))
  R=r+z;
  plot(R(1,:),R(2,:), 'g');
  hold on;  
end


%alpha1=a1(1:4:end);
%alpha2=a2(1:4:end);
%
%ar = (pi*R+L)/(2*R) - alpha1./R;
%
%x=(R + g_a.*cos(g_f.*ar)).*cos(ar)
%y=(R + g_a.*cos(g_f.*ar)).*sin(ar)
%
%r11=(2*g_a*g_f*K*K).*sin(g_f.*ar).*sin(ar)-(g_a*K*K.*cos(g_f.*ar).*(g_f*g_f+1)+1).*cos(ar);
%r12=(-2*g_a*g_f*K*K).*sin(g_f.*ar).*cos(ar)-(g_a*K*K.*cos(g_f.*ar).*(g_f*g_f+1)+1).*sin(ar);
%
alpha2 = a2(3);

x=r(1,:)+alpha2*n(1,:);
y=r(2,:)+alpha2*n(2,:);

l2 = a.*a+b.*b;
dl = sqrt(l2.*l2.*l2);

dn(1,:)=l2.*(2*g_a*g_f*K*K.*sin(g_f*ar).*cos(ar)+(g_a*K*K*(1+g_f*g_f).*cos(g_f*ar)+K).*sin(ar));
dn(2,:)=l2.*(2*g_a*g_f*K*K.*sin(g_f*ar).*sin(ar)-(g_a*K*K*(1+g_f*g_f).*cos(g_f*ar)+K).*cos(ar));
dn(1,:)=dn(1,:)-(1+g_a*(1-g_f*g_f)*K.*cos(g_f*ar)).*(g_a*g_f*K*K.*sin(g_f*ar)).*((1+g_a*K.*cos(g_f.*ar)).*cos(ar)-g_a*g_f*K*sin(g_f.*ar).*sin(ar));
dn(2,:)=dn(2,:)-(1+g_a*(1-g_f*g_f)*K.*cos(g_f*ar)).*(g_a*g_f*K*K.*sin(g_f*ar)).*((1+g_a*K.*cos(g_f.*ar)).*sin(ar)+g_a*g_f*K*sin(g_f.*ar).*cos(ar));
dn=dn./dl;

dr(1,:)=a.*cos(ar)+b.*sin(ar);
dr(2,:)=a.*sin(ar)-b.*cos(ar);

R1=dr+alpha2*dn;
R2=n;

ip=1;

%quiver (x(1:ip:end),y(1:ip:end), v(1, 1:ip:end), v(2, 1:ip:end), 'k')

%quiver (x(1:ip:end),y(1:ip:end), dr(1, 1:ip:end), dr(2, 1:ip:end), 'm')
%quiver (x(1:ip:end),y(1:ip:end), dn(1, 1:ip:end), dn(2, 1:ip:end), 'c')

quiver (x(1:ip:end),y(1:ip:end), R1(1, 1:ip:end), R1(2, 1:ip:end), 'r')
quiver (x(1:ip:end),y(1:ip:end), R2(1, 1:ip:end), R2(2, 1:ip:end), 'b')

daspect([1 1 1])
