function [ oArgs ] = StiffnessMatrix(model, geom, N, M)
  l = geom(1);
	curvature = geom(2);
  gA=geom(3);
  gV=geom(4);
  
  hDown = model(1);
  hTop = model(2);
  rho = model(3);
  E=model(4);
  v=model(5);

  feCount = N*M;
  count = 8*feCount;
  hDelta = (hTop - hDown)/M;
  lDelta = l/N;
	
  HardMatrix = zeros(2*(N+1)*(M+1), 2*(N+1)*(M+1));

	for i=1:feCount
    alpha1start = (rem(i,N))*lDelta;
    alpha1end = alpha1start + lDelta;
    alpha2end = hTop-(fix((i-1)/N))*hDelta
    alpha2start = alpha2end - hDelta
    
		localMatrix = GetLocalStiffnessMatrix(l, curvature, gA, gV, E, v, alpha1start, alpha1end, alpha2start, alpha2end);
		HardMatrix = SumMatrix(HardMatrix, localMatrix, i, N, M);
	end
  
	oArgs = HardMatrix;
  
end

