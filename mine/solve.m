function [eigvec eigval x]=solve(iArgs, staticIndecies)
%[h,l,K,rho,E,v,N];
  h = iArgs(1);
	l = iArgs(2);
	K = iArgs(3);
  rho = iArgs(4);
	E = iArgs(5);
	v = iArgs(6);
	N = iArgs(7);
  
  [meshBegin, meshEnd] = GenerateMesh(N, l);
  
  x = [meshBegin meshEnd(N)];
  
  s = StiffnessMatrix(iArgs, meshBegin, meshEnd)
  m = rho*MassMatrix(iArgs, meshBegin, meshEnd)
  s = applyStaticBoundaryConditionsToMatrix(s, staticIndecies);
  m = applyStaticBoundaryConditionsToMatrix(m, staticIndecies);
  [vec, lam] = eig (s, m);
  [vec, lam] = sortResults(vec, lam);
  
  eigvec = extendresultWithStaticBoundaryConditions(vec, staticIndecies);
  eigval = lam;
end