function [ oArgs ] = GetLocalMassMatrix(E, v, K, hDelta, lDelta, feIndex, N, hTop)
alpha1start = (rem(feIndex,N)-1)*lDelta;
alpha1end = alpha1start + lDelta;
alpha2start = hTop-(fix(feIndex/N) + 1)*hDelta;
alpha2end = alpha2start + hDelta;

f = @(psi, teta) localMassMatrix(psi, teta, K, alpha1start, alpha1end, alpha2start, alpha2end);

oArgs=quadgch5nodes2dim(f);
end

