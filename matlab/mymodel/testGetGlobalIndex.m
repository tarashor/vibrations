clear;
clc;

N=9;
M=6;
feIndex = 54;

for i=1:8
    for j=1:8
      gI = getGlobalIndex(i, feIndex, N, M);
      gJ = getGlobalIndex(j, feIndex, N, M);
      printf ("[ %d,  %d ] => [ %d,  %d ]\n",   i, j, gI, gJ);
      
    end
  end
  
  
%getGlobalIndex(7, 54, N, M)