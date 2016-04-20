function [med indMed minDist] = medoid(X)
% [med indMed minDist] = medoid(X)
%
% Computes the medoid out of a set of vectors.
% The medoid itself is returned as "med" and its index as "indMed".
%
% The Medoid is defined the vector that has the least total
% distance "minDist" (sum of distances) from all the other vectors in the
% said set. So the difference from the set Mean is effectively
% that the medoid is necessarily a member of the set; also,
% the mean will coincide with the medoid, should the mean be
% a member of the set.
%
% The columns of the input matrix X are the vectors in question,
% so X is size JxN where J are the number of our vector variates
% and N and is the number of vectors.
%
% G.Sfikas 17 Fev 2009
N = size(X, 2);
distanceMatrix = zeros(N);

for i = 1:N
    distanceMatrix(:, i) = squareDist(X, X(:, i), eye(size(X, 1)));
end
[minDist indMed] = min(sum(distanceMatrix));
med = X(:, indMed);
return;