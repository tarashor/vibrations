clear;
clc;
a=1;
b=1;
h=0.0003;
rho=2700;
E=69*1000000000;
v=0.33;

D=E*h*h*h/(12*(1-v*v));

aN = 50;
bN = 50;
ad=0:a/aN:a;
bd=0:b/bN:b;
t=10;



for i=1:aN+1
  for j=1:bN+1
    w(i,j)=-rectangleplate(ad(i), bd(j), t, a, b, D, rho, h);
  endfor
endfor

