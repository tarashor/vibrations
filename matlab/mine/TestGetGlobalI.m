clc;
clear;
N=20;
K=3;

dimOneLayer = 3*2*(N+1);
dimAllLayers = (2*K+1)*2*(N+1);

c = [];

for k=1:K
  for i=1:dimOneLayer
    e = fix((i-1)/6);
    rem_i = rem(i-1,6);
    s=0;
    if (rem_i == 0)
      s=10;
    elseif (rem_i == 1)
      s=11;
    elseif (rem_i == 2)
      s=12;
    elseif (rem_i == 3)
      s=30;
    elseif (rem_i == 4)
      s=31;
    elseif (rem_i == 5)
      s=32;
    endif
    
    i_new = getGlobalI(i, k, K);
    
    c = [c; [k, s, e, i_new]];
    
  end
end

c = sortrows(c, 4);

[rows col] = size(c);

for j=1:rows
  printf("U(%d) %d%d => %d\n", c(j, 1), c(j, 2), c(j, 3), c(j, 4));
end


getBoundaryConditionIndiciesForLayeredMatrix(N, K)
