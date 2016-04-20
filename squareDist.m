function res = squareDist(X, m, correl)
% Computes mahalanobis distance of each data (column) in X, to mean m,
% under correlation matrix correl. That is
%       x'* correl * x.
%
% See also mahalanobis.
%
% G.Sfikas 13 Nov 2007.
N = size(X, 2);
bigM = m * ones(1, N);
res = sum((X-bigM)' * (correl+eps*eye(size(X,1))) .* (X-bigM)', 2);
return;