function mod=applyStaticBoundaryConditionsToMatrix(M, N)
  indeciesToDelete = getBoundaryConditionIndicies(N);
  indecies=1:6*(N+1);
  indecies(indeciesToDelete)=[];
  mod=M(indecies, indecies);
end