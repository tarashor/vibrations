function [ Matrix ] = SumMatrix( globalMatrix, localMatrix, feIndex, N, M)

  for i=1:8
    for j=1:8
      gI = getGlobalIndex(i, feIndex, N, M);
      gJ = getGlobalIndex(j, feIndex, N, M);

      globalMatrix(gI, gJ) += localMatrix(i, j);
    end
  end
  Matrix=globalMatrix;
end

