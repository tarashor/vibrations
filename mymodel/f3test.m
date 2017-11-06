clear;
clc;
g_a = 0.03;
g_f = 20;
K=0.8;

R = 1/K;

L = 2;

p = 50;


delta = L/p;

a1 = 0:delta:L;
a2 = 0;

ar = pi/2+(L-2.*a1)/(2*R);

x=(R + g_a.*cos(g_f.*ar)).*cos(ar)
y=(R + g_a.*cos(g_f.*ar)).*sin(ar)

%ar_grad = ar*180/pi
figure
plot(x,y, 'g');

hold on;
alpha1=a1(1:4:end);
alpha2=a2(1:4:end);

ar = (pi*R+L)/(2*R) - alpha1./R;

x=(R + g_a.*cos(g_f.*ar)).*cos(ar)
y=(R + g_a.*cos(g_f.*ar)).*sin(ar)

%r11=(2*g_a*g_f*K*K).*sin(g_f.*ar).*sin(ar)-(g_a*K*K.*cos(g_f.*ar).*(g_f*g_f+1)+1).*cos(ar);
%r12=(-2*g_a*g_f*K*K).*sin(g_f.*ar).*cos(ar)-(g_a*K*K.*cos(g_f.*ar).*(g_f*g_f+1)+1).*sin(ar);

r11=(w.*sin(ar)+z.*cos(ar))
r12=(-w.*cos(ar)+z.*sin(ar))

quiver (x,y, r11, r12, 'r')

daspect([1 1 1])
