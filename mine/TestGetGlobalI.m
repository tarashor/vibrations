clc;
clear;
N=2;
K=3;

dimOneLayer = 3*2*(N+1);
dimAllLayers = (2*K+1)*2*(N+1);


for k=1:K
  for i=1:dimOneLayer
    e = fix((i-1)/6);
    rem_i = rem(i-1,6);
    s="";
    if (rem_i == 0)
      s="10";
    elseif (rem_i == 1)
      s="11";
    elseif (rem_i == 2)
      s="12";
    elseif (rem_i == 3)
      s="30";
    elseif (rem_i == 4)
      s="31";
    elseif (rem_i == 5)
      s="32";
    endif
    
    i_new = getGlobalI(i, k, K);
    printf("U(%d) %s%d => %d\n", k, s, e, i_new);
  end
end