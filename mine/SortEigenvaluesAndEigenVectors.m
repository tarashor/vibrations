function [v l]=SortEigenvaluesAndEigenVectors(vin, lin)
  lin = diag(lin);
  l=sort(lin);
  n=length(l);
  for i=1:n
    orig(i) = find(lin == l(i));
  end
  v=vin(:, orig);
end