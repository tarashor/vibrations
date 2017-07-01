function [ oArgs ] = StiffnessMatrix(model, geom, N, M)

  l = geom(1);
	curvature = geom(2);

  feCount = N*M;
  count = 8*feCount;
  hDelta = (model(2) - model(1))/M;
  lDelta = l/N;
  E=model(4);
  v=model(5);
	
  HardMatrix = zeros(count, count);

	for i=1:feCount
		localMatrix = GetLocalStiffnessMatrix(E, v, curvature, hDelta, lDelta, i, N, model(2));
		HardMatrix = SumMatrix(HardMatrix, localMatrix, i, N, M);
	end
	oArgs = HardMatrix;
end

