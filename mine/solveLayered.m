function [eigvec eigval x]=solveLayered(geom, layers, N, staticIndecies)  
  l = geom(1);
	curvature = geom(2);
  
  [meshBegin, meshEnd] = GenerateMesh(N, l);
  x = [meshBegin meshEnd(N)];

  layersCount = length(layers)
  
  dim = 2*(N+1)*(2*layersCount+1);
  SMatrix = zeros(dim,dim);
  MMatrix = zeros(dim,dim);

  for k=1:layersCount
    s = StiffnessMatrix(layers{k}, curvature, N, meshBegin, meshEnd);
    m = MassMatrix(layers{k}, N, meshBegin, meshEnd);
    
    SMatrix = SumLayerMatrix(SMatrix, s, k, layersCount);
    MMatrix = SumLayerMatrix(MMatrix, m, k, layersCount);
  end


  SMatrix = applyStaticBoundaryConditionsToMatrix(SMatrix, staticIndecies);
  MMatrix = applyStaticBoundaryConditionsToMatrix(MMatrix, staticIndecies);

  [vec, lam] = eig (SMatrix,MMatrix);
  [vec, lam] = SortEigenvaluesAndEigenVectors(vec, lam);

  eigvec = extendresultWithStaticBoundaryConditions(vec, staticIndecies);
  eigval = lam;
end