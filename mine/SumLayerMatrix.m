function [Matrix] = SumLayerMatrix(GlobalMatrix, LayerMatrix, k, layersCount)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
  [rows columns] = size(LayerMatrix);
  for i=1:rows
      i_new=getGlobalI(i, k, layersCount);
      for j=1:columns
        j_new=getGlobalI(j, k, layersCount);
        temp = GlobalMatrix(i_new, j_new);
        GlobalMatrix(i_new, j_new) = temp + LayerMatrix(i, j);
      end
  end
  Matrix=GlobalMatrix;
end