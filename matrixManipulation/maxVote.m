function y = maxVote(A)
% Use maximum-vote filter on K-level matrix 'A'.
% The filter used is 3x3.
%
% G. Sfikas 1/6/2006
%
K = max(max(A));
for j = -1:1
    for k = -1:1
        Neighbours((j+1)*3 + k+2,:) = neighbour(A, 1:(size(A,1)*size(A,2)), j, k);
    end
end
[junk winners] = max(hist(Neighbours, 0:K));
B = reshape(winners - 1, size(A));
y = A;
y(2:(size(A,1)-1), 2:(size(A,2)-1)) = B(2:(size(A,1)-1), 2:(size(A,2)-1));
return;

function y = neighbour(A, index, j, k)
% Find neighbour value of 'index' pixel in matrix A, which
% is +j, +k away from 'index'.
% ie neighbour(A, [4 100 55], -1, +1)
% will give NE neighbours of pixels with index 4, 100, 55.
% G. Sfikas 30/5/2006
B = zeros(size(A)+2);
B((2-j):(size(A,1)-j+1), (2-k):(size(A,2)-k+1)) = A;
C = B(2:(size(B,1)-1), 2:(size(B,2)-1));
y = C(index);
return;
