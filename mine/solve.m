function [eigvec eigval x]=solve(iArgs, staticIndecies)
%[h0,h1,l,K,rho,E,v,N];
  h0 = iArgs(1);
  h1 = iArgs(2);
	l = iArgs(3);
	K = iArgs(4);
  rho = iArgs(5);
	E = iArgs(6);
	v = iArgs(7);
	N = iArgs(8);
  
  [meshBegin, meshEnd] = GenerateMesh(N, l);
  
  x = [meshBegin meshEnd(N)];
  s = StiffnessMatrix(iArgs, meshBegin, meshEnd);
  m = rho*MassMatrix(iArgs, meshBegin, meshEnd);
  s = applyStaticBoundaryConditionsToMatrix(s, staticIndecies);
  m = applyStaticBoundaryConditionsToMatrix(m, staticIndecies);
  [vec, lam] = eig (s, m);
  [vec, lam] = sortResults(vec, lam);
  
  eigvec = extendresultWithStaticBoundaryConditions(vec, staticIndecies);
  eigval = lam;
end