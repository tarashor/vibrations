function [ oArgs ] = MassMatrix(model, geom, N, M)

  rho = model(3);
  l = geom(1);
	curvature = geom(2);

  feCount = N*M;
  count = 8*feCount;
  hDelta = (model(2) - model(1))/M;
  lDelta = l/N;
  E=model(4);
  v=model(5);
	
  MassMatrix = zeros(2*(N+1)*(M+1), 2*(N+1)*(M+1));

	for i=1:feCount
		localMatrix = GetLocalMassMatrix(E, v, curvature, hDelta, lDelta, i, N, model(2));
		MassMatrix = SumMatrix(MassMatrix, localMatrix, i, N, M);
	end
	oArgs = MassMatrix;
  
end
