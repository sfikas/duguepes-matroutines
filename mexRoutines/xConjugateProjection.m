function res = xConjugateProjection(x, weights)
% res = conjugateProjection(x, weights)
%
% Project Kx1 vector weights .* x to hyperplane 
% defined by weights'*y = 1, subject 
% to the constraint sum(y) = 1, y >= 0.
% The objective function is supposed to be of the form
%       x'Ax + bx + const.
% 
% Notes: 1. If weights = ones(1,K), then the result should
%           coincide to K.Blekas' BIDProjection(x).
% Notes: 2. In the CVPR08 submission, the objective function was
%           in fact of the form x'Ax + bx + clogx + const.
%           In BLEKAS05, the objective function function was
%           x'x + bx + clogx + const.
% G.Sfikas 21 April 2008

a = x(:);
K = numel(a);
weights = weights(:);
% Compute a useful function of 'weights'
g = weights.^(-1);
g = -g / sum(g);
% Initial projection to weights'y = 1.
y = a - g + g*sum(a);
yInit = y;
activeSet = zeros(K, 1);
for i = 1:K
    if(sum(y < 0) == 0)
        res = y;
        return;
    end
    activeSet = activeSet | (y < 0);
    lambda = zeros(K, 1);
    %%% Iterate between recalculating Lagrange operators & proj to sum(y) = 1
    for j = 1:K
        if activeSet(j) == 0
            lambda(j) = 0;
        else
            lambda(j) = -(yInit(j) + (g(activeSet)'*yInit(activeSet)) / sum(g(logical(1-activeSet))));
        end
    end
    for j = 1:K
        if activeSet(j) == 0
            y(j) = a(j) - g(j) + g(j)*sum(a) + g(j)*sum(lambda) + lambda(j);
        else
            y(j) = 0;
        end
    end
end
fprintf('conjugateProjection: Algorithm failed after K = %d iterations!', K);
return;