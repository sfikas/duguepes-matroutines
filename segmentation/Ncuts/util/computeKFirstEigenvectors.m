function [X,lambda,timing,D] = computeKFirstEigenvectors(W,nbEigenvectors,isncut,D);
% computes first k (normlalized-cuts) eigenvectors of W
% input:
% W: nxn affinity matrix
% nbEigenvectors: # eigenvectors requested (k)
% isncut: 0 for eigenvectors, 1 for normlalized-cuts eigenvectors
% D: optional parameter, can be set to diagonal of degree matrix of W; it
% is computed when not specified
% X: nxk eigenvectors
% lambda: kx1 eigenvalues
% timing: timing information
% typical calling sequence:
% [X,lambda] = computeKFirstEigenvectors(W,nbEigenvectors,1);
% Timothee Cour, 29-Aug-2006 07:49:15


if nargin < 3
    isncut = 1;
end

n = size(W,1);

[options,nbEigenvectors] = getDefaultOptionsEigs(n,nbEigenvectors);

if isncut
    if nargin>=4
        Dinvsqrt = 1./sqrt(D+eps);
        W=normalizeW_D(W,Dinvsqrt);
    else
        [W,Dinvsqrt,D]=normalizeW_D(W,[],0);
    end
end

if issparse(W)
    [result,timing] = eigs_optimized(@mex_w_times_x_symmetric_tril,[],nbEigenvectors,options,tril(W));
else
    [result,timing] = eigs_optimized(W,[],nbEigenvectors,options);
end
X = result.X;
if isncut
    X = spdiags(Dinvsqrt,0,n,n) * X;
end

lambda = result.lambda;
mex_normalizeColumns(X);
