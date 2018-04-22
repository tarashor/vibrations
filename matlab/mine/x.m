function [x] = x(resultVector, xVector, alfa, h)
%X Summary of this function goes here
%   Detailed explanation goes here
u30 = U30(resultVector, xVector,alfa(1)) * p0(alfa(3),h);
u31 = U31(resultVector, xVector,alfa(1)) * p1(alfa(3),h);
u32 = U32(resultVector, xVector,alfa(1)) * p2(alfa(3),h);
x = u30+u31+u32;  
%alfa(1)        
%k = FindCurIndex(xVector, alfa(1))
%x=0;
end

function [k]=FindCurIndex(xVector, alfa)
N=length(xVector);
L=xVector(N);
res = 1;
if ((alfa <= L) && (alfa >= 0))
   res = 2;
   while ((res < N) && (xVector(res) < alfa))      
        res = res + 1;
   end
end            
k=res;
end


function [u] = U30(resultVector, xVector, alfa)
    k = FindCurIndex(xVector, alfa);
    u = resultVector(6 * (k-1) - 2) * (xVector(k) - alfa) / (xVector(k) - xVector(k-1)) + resultVector(6 * (k-1) + 4) * (alfa - xVector(k-1)) / (xVector(k) - xVector(k-1));
end


function [u] = U31(resultVector, xVector, alfa)
    k = FindCurIndex(xVector, alfa);
    u = resultVector(6 * (k-1) - 1) * (xVector(k) - alfa) / (xVector(k) - xVector(k-1)) + resultVector(6 * (k-1) + 5) * (alfa - xVector(k-1)) / (xVector(k) - xVector(k-1));
end


function [u] = U32(resultVector, xVector, alfa)
    k = FindCurIndex(xVector, alfa);
    u = resultVector(6 * (k-1)) * (xVector(k) - alfa) / (xVector(k) - xVector(k-1)) + resultVector(6 * (k-1) + 6) * (alfa - xVector(k-1)) / (xVector(k) - xVector(k-1));
end



function [p]=p0(alfa3, h)
   p = 0.5 - alfa3 / h;
end

function [p]=p1(alfa3, h)
   p = 0.5 + alfa3 / h;
end
function [p]=p2(alfa3, h)
   p = 1 - (4 * alfa3 * alfa3) / (h * h);
end