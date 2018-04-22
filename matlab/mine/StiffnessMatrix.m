function [ oArgs ] = StiffnessMatrix(model, curvature, N, meshBegin, meshEnd)
	A = MatrixA(model);
	D = MatrixD(model, curvature);

	M = D' * A * D;

	count = 6 * (N + 1);

	HardMatrix = zeros(count, count);

	for i=1:N
		localMatrix = GetLocalMatrix(M, meshBegin(i), meshEnd(i));
		HardMatrix = SumMatrix(HardMatrix, localMatrix, i);
	end
	oArgs = HardMatrix;
end

