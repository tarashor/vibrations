function ind=getBoundaryConditionIndiciesForLayeredMatrix(N, layersCount)
  ind = [1,2*layersCount+2,(2*layersCount+1)*2*N+1, (2*layersCount+1)*(2*N+1)+1];
end