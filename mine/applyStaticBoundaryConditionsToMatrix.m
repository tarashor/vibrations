function mod=applyStaticBoundaryConditionsToMatrix(M, staticIndecies)
  indeciesToDelete = staticIndecies;
  indecies=1:length(M);
  indecies(indeciesToDelete)=[];
  mod=M(indecies, indecies);
end