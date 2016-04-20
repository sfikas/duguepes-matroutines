function res = gaussianValue(X, m, covar)
% Computes probability of dataset X (each _column_ is one datum),
% given a gaussian pdf, Normal(m, covar).
%
%   Example:
%           gaussianValue([4 3 2; 2 3 4], [3 3]', eye(2))
%       will return [0.0585 0.1592 0.0585]'.
%
% G.Sfikas 7 feb 2007.
N = size(X, 2);
d = size(X, 1);
bigM = m * ones(1, N);
mah = sum((X-bigM)' * (inv(covar)+eps*eye(size(X,1))) .* (X-bigM)', 2);
res = (2*pi)^(-0.5*d) * det(covar)^(-0.5) * exp(-0.5*mah);
return;