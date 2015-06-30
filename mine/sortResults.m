function [v l]=sortResults(vin, lin)
  lin = diag(lin);
  l=sort(lin);
  n=length(l);
  v=[];
  for i=1:n
    orig(i) = find(lin == l(i));
  end
  v=vin(:, orig);
end