function [ oArgs ] = MassMatrix( model, N, meshBegin, meshEnd )
  h0 = model(1);
  h1 = model(2);
  rho = model(3);
	E = model(4);
	v = model(5);
  h=h1-h0;

	count = 6 * (N + 1);

	MassMatrix = zeros(count, count);

  M = rho*baseAlfa3Functions(h);
	for i=1:N
		localMatrix = GetLocalMassMatrix(M, meshBegin(i), meshEnd(i));
		MassMatrix = SumMatrix(MassMatrix, localMatrix, i);
	end
	oArgs = MassMatrix;
end
