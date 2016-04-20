function [centroids weights idx covars] = deterministicKmeans(X, K)
% deterministicKmeans(X, K)
% Computes K centroids for data in X (each column is one datum), using
% the k-means algorithm. Uses a deterministic way, depending only on X,
% to initialize the k-means algorithm.
%
% Important notes:
%       The kmeans function in MATLAB2006b assumes that input is
%       a matrix where each _row_ is one datum.
%
% Arguments:
% X             Data. Each data is one column. (dxN)
% K             Number of desired centroids. (scalar)
%
% centroids     Resulting centroids. Each centroid is one column. (dxK)
% weights       Percentages showing how many data 'belong' to each
%                   centroid.
% idx           Segmentation (assignment of a cluster to each datum)
% covars        Covariance matrices for each cluster.
% 
%
% See also:
%       kmeans, multistartKmeans
% G.Sfikas 24 Apr 2007, 31 Jan 2008 (segmentation)
%          28 Apr 2008 (check for condition of output covariances)
%
weights = inv(K) * ones(1, K);
[d N] = size(X);
covars = zeros(d, d, K);
allDataMean = mean(X, 2);
allDataCov = cov(X');
centroids = zeros(d, K);
deviations = sqrt(diag(allDataCov));
if K == 1
    centroids(:, 1) = allDataMean;
else
    for i = 1:d
        centroids(i, :) = -deviations(i):(deviations(i)*2)/(K-1):deviations(i);
        centroids(i, :) = centroids(i, :) + allDataMean(d);
    end
    [idx, centroids] = kmeans(X', K, 'start', centroids', ...
        'EmptyAction', 'singleton', 'Maxiter', 200);
    centroids = centroids';
end

for j = 1:K
    weights(j) = sum(idx == j);
    covars(:,:,j) = cov(X(:,idx == j)') + eps*eye(d);
    if(rcond(covars(:,:,j)) <= 2*eps)
        covars(:,:,j) = eye(d);
    end
end
weights = weights / sum(weights);
return;