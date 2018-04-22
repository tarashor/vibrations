function [eigvec eigval]=solve(geom, model, N, M, staticIndecies)  

  SMatrix = StiffnessMatrix(model, geom, N, M);
  MMatrix = MassMatrix(model, geom, N, M);

  SMatrix = applyStaticBoundaryConditionsToMatrix(SMatrix, staticIndecies);
  MMatrix = applyStaticBoundaryConditionsToMatrix(MMatrix, staticIndecies);

  [vec, lam] = eig (SMatrix,MMatrix);
  
  [vec, lam] = SortEigenvaluesAndEigenVectors(vec, lam);

  eigvec = extendresultWithStaticBoundaryConditions(vec, staticIndecies);
  eigval = lam;
end