function [ meshBegin, meshEnd ] = GenerateMesh( N, l )
h=l/N;
meshBegin=[0];
meshEnd=[h];
for i=2:N
    meshBegin(i) = meshEnd(i-1);
    meshEnd(i)= meshBegin(i) + h;
end

end

