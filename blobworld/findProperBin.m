function y = findProperBin(A)
% Find proper bins for a 3xN matrix containing [0..255] values.
% L* is divided in 5 bins,
% a* and b* in 10 bins each.
% G.Sfikas 30/5/2006
A = double(A);
A(1,:) = A(1,:) / double(256/5);
A(2,:) = A(2,:) / double(256/10);
A(3,:) = A(3,:) / double(256/10);
A = floor(A);
y = floor([100 10 1] * double(A));
if (max(y) > 499)
    disp('Oh shit');
end
return;