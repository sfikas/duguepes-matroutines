function res = getBlobWorldBinCorelations(colorWeight, textureWeight)
%
% This will compute the bin corelation matrix 'A'. (p.1032 @ 2002paper)
% Same-numbered bins _are_ related (degree 1)
% Adjacent bins are somewhat related (degree 0.5)
%
% Additionally a weight can be set to count color or texture more than
% the other.
% Note that the resulting matrix is _not_ positive definite.
%
% G.Sfikas   15/6/2006
% Update/fix 29/6/2006
%
if nargin == 0
    colorWeight = 1;
    textureWeight = 1;
elseif nargin == 1
    textureWeight = 1;
end    
A = double(zeros(500, 500));
for i = 1:500
    for j = (i+1):500
        iH = floor(i/100);
        iT = floor((i - iH*100)/10);
        iU = floor(i - iH*100 - iT*10);
        %
        jH = floor(j/100);
        jT = floor((j - jH*100)/10);
        jU = floor(j - jH*100 - jT*10);
        if abs(iH - jH) + abs(iT - jT) + abs(iU - jU) == 1
            A(i,j) = 0.5;
        end
    end
end
A = A + A' + eye(500);
res = colorWeight * A;
res(501:502, 501:502) = textureWeight * eye(2);
return;