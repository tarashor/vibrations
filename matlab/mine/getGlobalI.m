function [i_new] = getGlobalI(i, k, layersCount)
      r_i=rem((i-1), 6);
      delta_i=r_i+1;
      if ((r_i<=5)&&(r_i>=3))
        delta_i = (2*layersCount+1)+r_i-2;
      endif
      l_i=fix((i-1)./6)+1;
      i_new=(2*layersCount+1)*2*(l_i-1)+(k-1)*2;
      i_new=i_new+delta_i;
end