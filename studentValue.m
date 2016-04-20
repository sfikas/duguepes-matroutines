function res = studentValue(X, m, covar, v)
% Computes probability of dataset X (each _column_ is one datum),
% given a student-t pdf, Student(m, covar, v).
%
%   Example:
%           studentValue([4 3 2; 2 3 4], [3 3]', eye(2), 2)
%       will return [0.0398 0.1592 0.0398]'.
%
% G.Sfikas 7 feb 2007.
N = size(X, 2);
d = size(X, 1);
bigM = m * ones(1, N);
mah = sum((X-bigM)' * (inv(covar)+eps*eye(size(X,1))) .* (X-bigM)', 2);
res = gamma(0.5*(v + d)) * gamma(0.5*v)^-1 * (pi*v)^-(0.5*d) * ...
    det(covar)^-0.5 * (1 + mah/v).^(-0.5*(v + d));
return;