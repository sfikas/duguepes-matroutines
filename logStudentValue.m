function res = logStudentValue(X, m, covar, v)
% Computes log-probability of dataset X (each _column_ is one datum),
% given a student-t pdf, Student(m, covar, v).
%
% G.Sfikas 5 feb 2008. 
% Update log(1 + .) -> log1p 19 Feb 2008
N = size(X, 2);
d = size(X, 1);
bigM = m * ones(1, N);
mah = sum((X-bigM)' * (inv(covar)+eps*eye(size(X,1))) .* (X-bigM)', 2);
res = gammaln(0.5*(v + d)) - gammaln(0.5*v) -(0.5*d)*log(pi*v)  ...
    -0.5*logdet(covar) -0.5*(v + d).*log1p(mah/v);
return;