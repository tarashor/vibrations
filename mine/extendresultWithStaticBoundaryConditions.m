function [ext]=extendresultWithStaticBoundaryConditions(res, N)
  indeciesToExtend  = getBoundaryConditionIndicies(N);
  extLength = length(res)+length(indeciesToExtend);
  ext=zeros(extLength,1);
  z=0;
  for i=1:extLength
    if (length(find(indeciesToExtend==i))==0)
      ext(i)=res(i-z);
    else
      z = z+1;
      
    end
  end
end