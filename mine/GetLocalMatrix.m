function [ oArgs ] = GetLocalMatrix(M, b, e)
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
    m = 1;
    for j=1:12
        if (rem(j,2)==1)
            C(j, m) = (1 - eta) / 2;
            C(j, m + 6) = (1 + eta) / 2;
        else
            C(j, m) = 1 / (b - e);
            C(j, m + 6) = 1 / (e - b);
            m = m + 1;
        end
    end
    
    T = C'*M*C;
    
    localMatrix = localMatrix + gw(i)*T;
    
end
localMatrix = localMatrix*((e - b) / 2.0);
oArgs=localMatrix;
end

