function [ oArgs ] = StiffnessMatrix( iArgs, meshBegin, meshEnd )
	N = iArgs(8);

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

