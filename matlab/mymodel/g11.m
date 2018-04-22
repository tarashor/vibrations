function [ oArgs ] = g11(l, K, gA, gV, alpha1, alpha2)

q=1+K*alpha2;
a=(pi+K*l)/2-K*alpha1;
w = q+gA*K*cos(gV*a);
z = gA*gV*K*sin(gV*a);

oArgs = 1/(w*w+z*z);

end
