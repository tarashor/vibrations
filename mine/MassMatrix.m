function [ oArgs ] = MassMatrix( iArgs, meshBegin, meshEnd )
	h = iArgs(1);
	l = iArgs(2);
	K = iArgs(3);
	E = iArgs(4);
	v = iArgs(5);
	N = iArgs(6);

	count = 6 * (N + 1);

	MassMatrix = zeros(count, count);

  M = baseAlfa3Functions(h);
	for i=1:N
		localMatrix = GetLocalMassMatrix(M, meshBegin(i), meshEnd(i));
		MassMatrix = SumMatrix(MassMatrix, localMatrix, i);
	end
	oArgs = MassMatrix;
end
