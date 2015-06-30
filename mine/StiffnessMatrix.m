function [ oArgs ] = StiffnessMatrix( iArgs, meshBegin, meshEnd )
  h = iArgs(1);
	l = iArgs(2);
	K = iArgs(3);
  rho = iArgs(4);
	E = iArgs(5);
	v = iArgs(6);
	N = iArgs(7);

	A = MatrixA(iArgs);
	D = MatrixD(iArgs);

	M = D' * A * D;

	count = 6 * (N + 1);

	HardMatrix = zeros(count, count);

	for i=1:N
		localMatrix = GetLocalMatrix(M, meshBegin(i), meshEnd(i));
		HardMatrix = SumMatrix(HardMatrix, localMatrix, i);
	end
	oArgs = HardMatrix;
end

