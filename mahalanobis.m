function res = mahalanobis(X, m, covar)
% Computes mahalanobis distance of each data (column) in X, to mean m,
% under correlation matrix covar. That is
%       x'* inv(covar)* x.
%
% G.Sfikas 7 feb 2007.
N = size(X, 2);
bigM = m * ones(1, N);
res = sum((X-bigM)' * (inv(covar)+eps*eye(size(X,1))) .* (X-bigM)', 2);
return;