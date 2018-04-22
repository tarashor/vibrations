function [ oArgs ] = GetLocalMassMatrix(M, b, e)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
            

  gn(1) = -sqrt(3/5);
  gn(2) = 0;
  gn(3) = sqrt(3/5);

  gw(1) = 5/9;
  gw(2) = 8/9;
  gw(3) = 5/9;

  gorder = 3;

  localMatrix = zeros(12,12);

  for i=1:gorder
    eta = gn(i);
    for j=1:6
      C(j, j) = (1 - eta) / 2;
      C(j, j + 6) = (1 + eta) / 2;
    end
    
    
    T = C'*M*C;
    
    localMatrix = localMatrix + gw(i)*T;
  end
  localMatrix = localMatrix*((e - b) / 2.0);
  oArgs=localMatrix;
end

