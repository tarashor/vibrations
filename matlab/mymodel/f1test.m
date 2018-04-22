clear;
clc;

g_a = 0.03;
g_f = 20;

R = 1/0.8;

L = 2;
p = 50;



%L = R*2*start_corner;

delta = L/p;
figure
hold on;
for a2=-0.025:0.05:0.025

  a1 = 0:delta:L;

  ar = (pi*R+L)/(2*R) - a1/R;

  x=(R + a2).*cos(ar);
  y=(R + a2).*sin(ar);

  %ar_grad = ar*180/pi
  
  plot(x,y, 'g');
  alpha1=a1(1:4:end);
  alpha2=a2(1:4:end);

  q=1+alpha2./R;
  ar2 = (alpha1 - L/2)./R


  x=q.*sin(ar2).*R
  y=q.*cos(ar2).*R

  r11=cos(ar2)
  r12=-sin(ar2)

  r21=sin(ar2)
  r22=cos(ar2)

  quiver (x,y, r11, r12, 'r')
  quiver (x,y, r21,r22,'b')
 end
 
 daspect([1 1 1])
 hold off;