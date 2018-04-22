function [oArgs] = quadgch5nodes2dim(f)
N=5;
w=[0.23692689, 0.47862867, 0.56888889, 0.47862867, 0.23692689];
x=[-0.90617985, -0.53846931, 0, 0.53846931, 0.90617985];


res = w(1) * w(1) * f(x(1), x(1));

for i=1:N
  for j=1:N
    if (i != 1 || j != 1)
      res += w(i) * w(j) * f(x(i), x(j));
    end
  end
end

oArgs = res;

end