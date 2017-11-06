clear;
clc;
g_a = 0.03;
g_f = 4;

R = 1/0.8;

L = 2;

p = 50;


delta = L/p;

a1 = 0:delta:L;
a2 = 0.5;

ar = pi/2+(L-2.*a1)/(2*R);

x=R.*cos(ar)-a2.*sin(ar);
y=R.*sin(ar)+a2.*cos(ar);

%ar_grad = ar*180/pi
figure
plot(x,y, 'g');
daspect([1 1 1])


