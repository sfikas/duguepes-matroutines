function [m, covar, v, w, z, tt0_100] = studentMixEmFit(X, nKernels)
% Maximum-likelihood fit of data 'X' on a student-t mixture model,
% made up by 'nKernels' number of kernels. All variates data vectors
% (X's columns) are assumed to be observed.
%
% Example:
%       studentMixEmFit(randn(2, 100));
%   will fit 100 random bivariate data on a 4-kernel student-t mixture.
%       studentMixEmFit(X, 11);
%   will fit data from columns of matrix 'X' on a 11-kernel student-t
%   mixture.
%
% Arguments:
% X         -   Input. Each column is one datum.
% nKernels  -   Number of kernels on mixture model.
% m         -   EM means; i-th column represents mean of the i-th kernel.
% covar     -   EM covariances; i-th matrix represents covariance of the
%                   i-th kernel. (ie (:,:,i) ).
% v         -   EM degrees of freedom.
% w         -   EM kernel weigths.
%
% G.Sfikas 29.11.2006
%
% Update 10 feb 2007 (print average logLikelihood)
% Update 26 Nov 2007
d = size(X, 1);
N = size(X, 2);
if nargin < 2
    nKernels = 4;
end
% allDataMean = mean(X, 2);
% allDataCov = cov(X');
% % Random initialization:
% % m = sqrtm(allDataCov) * randn(d, nKernels) + allDataMean * ones(1, nKernels);
% % covar = zeros(d, d, nKernels);
% % for i = 1:nKernels
% %     covar(:,:,i) = allDataCov / nKernels^2;
% % end
% % FIN - Random initialization
% % Deterministic initialization -- data independent
% % m = zeros(d, nKernels);
% % for i = 1:nKernels
% %     covar(:,:,i) = eye(d);
% % end
% % FIN - Deterministic initialization -- data independent
% % Deterministic initialization -- K-means
% % TODO: Replace this with a call to "deterministicKmeans"
% deviations = sqrt(diag(allDataCov));
% if nKernels == 1
%     m(:, 1) = allDataMean;
% else
%     for i = 1:d
%         m(i, :) = -deviations(i):(deviations(i)*2)/(nKernels-1):deviations(i);
%         m(i, :) = m(i, :) + allDataMean(d);        
%     end
%     [junk, m] = kmeans(X', nKernels, 'start', m', ...
%     'EmptyAction', 'singleton');
%     m = m';
% end
% for i = 1:nKernels
%     % Make each std.deviation equal to 1/nKernels of total std.deviation.
%     covar(:,:,i) = allDataCov / nKernels^2;
% end
% FIN - Deterministic initialization -- K-means
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% w = ones(1, nKernels) / nKernels;
%Ta varh tha mporousan na paroun times analogws me to plhthos twn data
% pou anhkoun se kathe pyrina, symfwna me thn arxikopoiish tou k-means.
[m w junk covar] = deterministicKmeans(X(:,1:N), nKernels);
disp('Warning: Initializing all kernels to v = 1.');
v = ones(1, nKernels) * 1;
clear allDataMean;
clear allDataCov;
% End Init
%
% Main EM loop
%
iteration = 0;
z = zeros(nKernels, N);
t = zeros(nKernels, N);
llChangeRatio = 0;
disp('EM iteration  Likelihood     Average Logl      LhoodIncrease');
disp('------------------------------------------------------------');
while 1
    % E-step
    for i = 1:nKernels
        z(i, :) = w(i) * computeProbability(X, m(:,i), covar(:,:,i), v(i));
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
    for i = 1:nKernels
        t(i, :) = (v(i) + d) ./ (v(i) + mahalanobis(X, m(:,i), covar(:,:,i)));
    end
    % M-step
    w = sum(z, 2)' / N;
    for i = 1:nKernels
        if sum(t(i, :).*z(i, :)) > 0 %To catch empty cluster
            m(:, i) = (X * (t(i, :) .* z(i, :))') / sum(t(i, :).*z(i, :));
        end
%       Normal EM covariance update:        
%         newCovar = inv(sum(z(i, :))) * ( ...        
%       [Meng,VanDyk 97] modification update for faster convergence:
        newCovar = inv(sum(z(i, :) .* t(i, :))) * ( ...
            (X - m(:, i)*ones(1,N))* ...
            sparse(1:N, 1:N, t(i, :).*z(i, :), N, N)* ...
            (X - m(:, i)*ones(1,N))' );
        newCovar = newCovar + eps*eye(d);
        if rcond(newCovar) > 1e3 * eps
            covar(:, :, i) = newCovar;
        end
        clear newCovar;
        vconstant = 1 + psi(0.5*(v(i) + d)) - log(0.5*(v(i) + d)) + ...
            inv(sum(z(i, :))) * (log(t(i, :)) - t(i, :)) * z(i, :)';
        v(i) = bisection(vconstant);
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
    tt0_100(j) = w(1) * computeProbability(tt0_100range(j), m(1), covar(:, :, 1), v(1)) + ...
        w(2) * computeProbability(tt0_100range(j), m(2), covar(:, :, 2), v(2));
    j = j + 1;
end
return;

function res = mahalanobis(X, m, covar)
N = size(X, 2);
bigM = m * ones(1, N);
res = sum((X-bigM)' * (inv(covar)+eps*eye(size(X,1))) .* (X-bigM)', 2);
return;

function res = computeProbability(X, m, covar, v)
% Compute probability of data in 'X', assuming it was generated by a
% student-t with parameters 'm', 'covar', 'v'.
d = size(X, 1);
t1 = gamma(0.5*(v + d)) / gamma(0.5*v);
t2 = det(covar)^(-0.5) * (pi*v)^(-0.5*d);
t3 = (1 + mahalanobis(X, m, covar)/v).^(-0.5*(v+d));
res = t1 * t2 * t3;
return;

function res = bisection(k)
%Init
leftbound = 10e-3;
rightbound = 10;
while 1
    if log(0.5*rightbound) - psi(0.5*rightbound) + k > 0
        rightbound = rightbound*2;
    else
        break;
    end
end
%Start
iter = 0;
while 1
    x = 0.5 * (leftbound + rightbound);
    y = log(0.5*x) - psi(0.5*x) + k;
    if abs(y) < 1e-5
        res = x;
        break;
    elseif y > 0
        leftbound = x;
    elseif y < 0
        rightbound = x;
    end
    iter = iter+1;
end %of while
return;

function res = computeLikelihood(X, m, covar, v, w)
nKernels = size(w, 2);
for i = 1:nKernels
        z(i, :) = w(i) * computeProbability(X, m(:,i), covar(:,:,i), v(i));
end
res = sum(log(sum(z, 1)));
return;