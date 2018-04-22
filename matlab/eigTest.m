clc;
clear;
a=[1 2; 3 4];
b = [5; 10];
[v l]=eig(a);

a*v(:,1)
l(1)*v(:,1)