function [eigvec eigval x]=solveOneLayer(geom, model, N, staticIndecies)
  l = geom(1);
	curvature = geom(2);
  
  [meshBegin, meshEnd] = GenerateMesh(N, l);
  x = [meshBegin meshEnd(N)];
  
  s = StiffnessMatrix(model, curvature, N, meshBegin, meshEnd);
  m = MassMatrix(model, N, meshBegin, meshEnd);
  
  s = applyStaticBoundaryConditionsToMatrix(s, staticIndecies);
  m = applyStaticBoundaryConditionsToMatrix(m, staticIndecies);
  
  [vec, lam] = eig (s, m);
  [vec, lam] = SortEigenvaluesAndEigenVectors(vec, lam);
  
  eigvec = extendresultWithStaticBoundaryConditions(vec, staticIndecies);
  eigval = lam;
end