function [x1 x2]=X(a1, a2, R, L)
x1=(a2+R)'*cos((pi*R+L-2*a1)./(2*R));
x2=(a2+R)'*sin((pi*R+L-2*a1)./(2*R));
end