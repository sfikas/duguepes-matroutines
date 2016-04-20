function res = logGaussianValue(X, m, covar)
% Computes log-probability of dataset X (each _column_ is one datum),
% given a gaussian pdf, Normal(m, covar).
%
% See also:
%       gaussianValue, studentValue
%
% G.Sfikas 4 feb 2008.
N = size(X, 2);
d = size(X, 1);
bigM = m * ones(1, N);
mah = sum((X-bigM)' * (inv(covar)+eps*eye(size(X,1))) .* (X-bigM)', 2);
res = -(0.5*d)*log(2*pi) -0.5*logdet(covar+eps) - 0.5*mah;
return;