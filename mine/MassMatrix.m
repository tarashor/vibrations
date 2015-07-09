function [ oArgs ] = MassMatrix( iArgs, meshBegin, meshEnd )
%[h0,h1,l,K,rho,E,v,N];
  h0 = iArgs(1);
  h1 = iArgs(2);
	l = iArgs(3);
	K = iArgs(4);
  rho = iArgs(5);
	E = iArgs(6);
	v = iArgs(7);
	N = iArgs(8);
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
