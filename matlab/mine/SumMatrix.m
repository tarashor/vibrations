function [ Matrix ] = SumMatrix( MMatrix, NMatrix, k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[rows columns] = size(NMatrix);
for i=1:rows
    for j=1:columns
        temp = MMatrix(6 * (k-1) + i, 6 * (k-1) + j);
        MMatrix(6 * (k-1) + i, 6 * (k-1) + j) = temp + NMatrix(i, j);
    end
end
Matrix=MMatrix;
end

