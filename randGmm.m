function X = randGmm(N, w, m, U, V)
% RANDGMM   Sample from mixture of multivariate normals.
% Use as X = randgmm(N, w, m, U)
% or     X = randgmm(N, w, m, [], V)
% 
% Returns a Nxd matrix, where each column is a sample from
% a mixture of multivariate normals. The parameters of the mixture are
% w     Kx1 kernel weight vector.
% m     dxK matrix of means (each mean is a column)
% U     dxdxK matrix of upper triangular Cholesky factors of the
%       covariance matrices.
% V     dxdxK matrix of covariance matrices.
%
% Requires "lightspeed" library.
% G.Sfikas 20 Mar 2007
%

d = size(m, 1);
K = size(m, 2);
X = [];
% Use a multinomial sample first, to decide how many data
% to 'assign' to each cluster.
% And make sure w is a column matrix.
if size(w, 2) > size(w, 1)
    w = w';
end
Np = sample_hist(w, N);
if nargin < 5
    for i = 1:K
        X = [X, randnorm(Np(i), m(:, i), U(:,:,i))];
    end
else
    for i = 1:K
        X = [X, randnorm(Np(i), m(:, i), [], V(:,:,i))];
    end
end

return;