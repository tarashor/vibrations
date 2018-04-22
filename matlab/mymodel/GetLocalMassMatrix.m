function [ oArgs ] = GetLocalMassMatrix(l, K, gA, gV, E, v, alpha1start, alpha1end, alpha2start, alpha2end)
  f = @(psi, teta) localMassMatrix(psi, teta, l, K, gA, gV, alpha1start, alpha1end, alpha2start, alpha2end);
  oArgs=quadgch5nodes2dim(f);
end

