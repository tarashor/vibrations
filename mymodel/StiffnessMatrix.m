function [ oArgs ] = StiffnessMatrix(model, geom, N, M)

  l = geom(1);
	curvature = geom(2);

  feCount = N*M;
  count = 8*feCount;
  hDelta = (model(2) - model(1))/M;
  lDelta = l/N;
	
  HardMatrix = zeros(count, count);

	for i=1:feCount
		localMatrix = GetLocalStiffnessMatrix(M, );
		HardMatrix = SumMatrix(HardMatrix, localMatrix, i);
	end
	oArgs = HardMatrix;
end

