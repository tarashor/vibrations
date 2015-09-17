function [ext]=extendresultWithStaticBoundaryConditions(v, staticIndecies)
  indeciesToExtend  = staticIndecies;
  extLength = size(v,1) + length(indeciesToExtend);
  ext=[];
  z=0;
  for i=1:extLength
    if (length(find(indeciesToExtend==i))==0)
      ext=[ext;v(i-z,:)];
    else
      z = z+1;
      ext=[ext; zeros(1,length(v))];
    end
  end
end