function [m, covar, w, z, tt0_100] = gaussianMixEmFit(X, nKernels)
% Maximum-likelihood fit of data 'X' on a gaussian mixture model,
% made up by 'nKernels' number of kernels. All variates data vectors
% (X's columns) are assumed to be observed.
%
% Example:
%       gaussianMixEmFit(randn(2, 100));
%   will fit 100 random bivariate data on a 4-kernel (default) gaussian
%   mixture.
%       gaussianMixEmFit(X, 11);
%   will fit data from columns of matrix 'X' on a 11-kernel gaussian
%   mixture.
%
% Arguments:
% X         -   Input. Each column is one datum.
% nKernels  -   Number of kernels on mixture model.
% m         -   EM means; i-th column represents mean of the i-th kernel.
% covar     -   EM covariances; i-th matrix represents covariance of the
%                   i-th kernel. (ie (:,:,i) ).
% w         -   EM kernel weigths.
%
% G.Sfikas 19.12.2006
%
% Update 10 feb 2007 (print average logLikelihood)
% Update 12 nov 2007 (added Z on output)
% Update 28 nov 2008 (use 'gaussianValue' sfikaslib function)
%
d = size(X, 1);
N = size(X, 2);
if nargin < 2
    nKernels = 4;
end
% Initialization
% allDataCov = cov(X');
% covar = zeros(d, d, nKernels);
% % Random initialization:
% % m = sqrtm(allDataCov) * randn(d, nKernels) + allDataMean * ones(1, nKernels);
% % covar = zeros(d, d, nKernels);
% % for i = 1:nKernels
% %     covar(:,:,i) = allDataCov / nKernels^2;
% % end
% % FIN - Random initialization
% % Deterministic initialization -- K-means
% [m w] = deterministicKmeans(X, nKernels);
% for i = 1:nKernels
%     % Make each std.deviation equal to 1/nKernels of total std.deviation.
%     covar(:,:,i) = allDataCov / nKernels^2 + eps*eye(d);
% end
% FIN - Deterministic initialization -- K-means
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m w junk covar] = deterministicKmeans(X(:,1:N), nKernels);
clear allDataMean;
clear allDataCov;


% End Init
%
% Main EM loop
%
iteration = 0;
llChangeRatio = 0;
disp('EM iteration  Likelihood     Average Logl      LhoodIncrease');
disp('------------------------------------------------------------');
while 1
    % E-step
    for i = 1:nKernels
        z(i, :) = w(i) * gaussianValue(X, m(:,i), covar(:,:,i));
    end
    % Before 'z' is correctly normalized, use it to compute the likelihood:
    likelihood = sum(log(sum(z, 1)+eps));
    z = z ./ (ones(nKernels, 1) * (sum(z, 1)+eps));
    if iteration ~= 0
        llChangeRatio = (likelihood - prev) / abs(prev);
    end
    disp(sprintf('%3d           %3.2f         %2.5f      %2.5f%%', ...
            iteration, likelihood, likelihood/N, llChangeRatio*100));
    if iteration > 1 && llChangeRatio < 1e-3
        break;
    end
    iteration = iteration + 1;
    % M-step
    w = sum(z, 2)' / N;
    for i = 1:nKernels
        if(sum(z(i, :)) > 0) %To catch empty clusters
            m(:, i) = (X * z(i, :)') / sum(z(i, :));
        end
        newCovar = inv(sum(z(i, :))) * ( ...
            (X - m(:, i)*ones(1,N))* ...
            sparse(1:N, 1:N, z(i,:), N, N) * ...
            (X - m(:, i)*ones(1,N))' );
        newCovar = newCovar + eps*eye(d);
        if rcond(newCovar) > 1e3 * eps
            covar(:, :, i) = newCovar;
        end
        clear newCovar;
    end
    prev = likelihood;
end
tt0_100range = 0:0.1:100;
tt0_100 = zeros(numel(tt0_100range), 1);
j = 1;
%v(2) = .1;
%w(2) = 2*w(2);
%w(1) = 1 - w(2);
for i = tt0_100range
    tt0_100(j) = w(1) *  mvnormpdf(tt0_100range(j), m(1), [], covar(:, :, 1)) + ...
        w(2) * mvnormpdf(tt0_100range(j), m(2), [], covar(:, :, 2));
    j = j + 1;
end
return;