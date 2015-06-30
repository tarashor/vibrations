function [ oArgs ] = MatrixA( iArgs )
%K Summary of this function goes here
%   Detailed explanation goes here

h = iArgs(1);
l = iArgs(2);
K = iArgs(3);
E = iArgs(4);
v = iArgs(5);

A(1, 1) = (E * h * (1 - v)) / (3 * (1 + v) * (1 - 2 * v));
A(1, 2) =(E * h * (1 - v)) / (6 * (1 + v) * (1 - 2 * v));
A(1, 3) = (E * h * (1 - v)) / (3 * (1 + v) * (1 - 2 * v));
A(1, 4) = (E * h * v) / (3 * (1 + v) * (1 - 2 * v));
A(1, 5) = (E * h * v) / (6 * (1 + v) * (1 - 2 * v));
A(1, 6) = (E * h * v) / (3 * (1 + v) * (1 - 2 * v));

A(2, 2) = (E * h * (1 - v)) / (3 * (1 + v) * (1 - 2 * v));

A(2, 3) = (E * h * (1 - v)) / (3 * (1 + v) * (1 - 2 * v));
A(3, 1) = A(2, 3);
A(3, 2) = A(2, 3);
A(4, 4) = A(2, 3);
A(4, 6) = A(2, 3);
A(5, 5) = A(2, 3);
A(5, 6) = A(2, 3);
A(6, 4) = A(2, 3);
A(6, 5) = A(2, 3);

A(2, 1) = (E * h * (1 - v)) / (6 * (1 + v) * (1 - 2 * v));
A(4, 5) = A(2, 1);
A(5, 4) = A(2, 1);

A(3, 3) = (8 * E * h * (1 - v)) / (5 * (1 + v) * (1 - 2 * v)); 
A(6, 6) = A(3, 3);


A(2, 5) = (E * h * v) / (3 * (1 + v) * (1 - 2 * v));
A(2, 6) = A(2, 5);
A(3, 4) = A(2, 5);
A(3, 5) = A(2, 5);
A(4, 1) = A(2, 5);
A(4, 3) = A(2, 5);
A(5, 2) = A(2, 5);
A(5, 3) = A(2, 5);
A(6, 1) = A(2, 5);
A(6, 2) = A(2, 5);
             
A(2, 4) = (E * h * v) / (6 * (1 + v) * (1 - 2 * v));
A(4, 2) = A(2, 4);
A(5, 1) = A(2, 4);

A(3, 6) = (8 * E * h * v) / (5 * (1 + v) * (1 - 2 * v));
A(6, 3) = A(3, 6);

A(7, 7) = (E * h) / (6 * (1 + v));
A(7, 9) = A(7, 7);
A(8, 8) = A(7, 7);
A(8, 9) = A(7, 7);
A(9, 7) = A(7, 7);
A(9, 8) = A(7, 7);

A(7, 8) = (E * h) / (12 * (1 + v));
A(8, 7) = A(7, 8);
A(9, 9) = (4 * E * h) / (15 * (1 + v));

oArgs = A;
end

