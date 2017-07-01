function [eigvec eigval x]=solve(geom, model, N, M, staticIndecies)  
  l = geom(1);
	curvature = geom(2);

  SMatrix = StiffnessMatrix(model, geom, N, M);
  MMatrix = MassMatrix(model, N, meshBegin, meshEnd);

  SMatrix = applyStaticBoundaryConditionsToMatrix(SMatrix, staticIndecies);
  MMatrix = applyStaticBoundaryConditionsToMatrix(MMatrix, staticIndecies);

  [vec, lam] = eig (SMatrix,MMatrix);
  
  [vec, lam] = SortEigenvaluesAndEigenVectors(vec, lam);

  eigvec = extendresultWithStaticBoundaryConditions(vec, staticIndecies);
  eigval = lam;
end