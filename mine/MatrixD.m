function [ oArgs ] = MatrixD( iArgs )
%K Summary of this function goes here
%   Detailed explanation goes here

h = iArgs(1);
l = iArgs(2);
K = iArgs(3);
E = iArgs(4);
v = iArgs(5);

hInv=1/h;

D(1, 2) = 1; 
D(2, 4) = 1;
D(3, 6) = 1;
D(7, 8) = 1;
D(8, 10) = 1;
D(9, 12) = 1;

D(1, 7) = K;
D(2, 9) = K;
D(3, 11) = K;
D(9, 5) = K;
D(6, 11) = 2 * K;

D(4, 7) = -hInv + K / 2;
D(4, 9) = hInv - K / 2;
D(7, 3) = hInv - K / 2;
D(8, 3) = hInv - K / 2;

D(4, 11) = 4 * hInv - 2 * K;
D(7, 5) = 4 * hInv - 2 * K;

D(5, 7) = -hInv - K / 2;
D(7, 1) = -hInv - K / 2;
D(8, 1) = -hInv - K / 2;
            
D(5, 9) = hInv + K / 2;

D(5, 11) = -4 * hInv - 2 * K;
D(8, 5) = -4 * hInv - 2 * K;

oArgs = D;
end

