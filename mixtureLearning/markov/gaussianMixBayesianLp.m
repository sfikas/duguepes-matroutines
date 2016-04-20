function [m, covar, w, z, u, a_ksi, beta, likelihood, w5, z5] = gaussianMixBayesianLp(X, K, imageSize, groundTruth, options, edgemap)
% [m, covar, w, z, u] = gaussianMixBayesianLp(X, maxSegments, imageSize);
% 
% Gaussian mixture model fit to data X. A _compound_ markov random field 
% is assumed  on the spatially-varying weights, with a Gaussian field clique 
% function. A first-order (4-neighbourhood) is considered on the MRF.
% The cMRF class variances (beta) are _directional_ and _class_ adaptive.
% We also assume the mean-field approximation
%           q(u, ksi) = q(u)q(ksi)
% The inference method is _Variational_ inference.
%
%       | 2 |
%     --/   \--
%     1   x   3     <-- Neighbourhood indexing
%     --\   /--
%       | 4 |
%
% Update history 
%       1. 13 May 2008
%       2. 25 Oct 2008 (warning message for adding noise on init)
if (exist('edgemap', 'var') == 0)
    disp('gaussianMixBayesianLp: Will work without an initial edgemap');
end
flagGridScanAlsoZ = 0;
flagGridScan = 0;
flagMultiRes = 0;
flagEdgeDetector = 0;
flagEdgeDetectorFIX = 0;
if(exist('options', 'var') == 1)
    for oCount = 1:numel(options)
        switch lower(options(oCount))
            case {'g'}
                flagGridScan = 1;
                disp('gaussianMixBayesianLp: Using grid-scan optimization on Pi');                
            case {'z'}
                flagGridScanAlsoZ = 1;
                disp('gaussianMixBayesianLp: Using grid-scan also on Z');                
            case {'m'}
                flagMultiRes = 1;
                disp('gaussianMixBayesianLp: Using multiresolution optimization');
            case {'e'}
                if flagEdgeDetector == 0
                    flagEdgeDetector = 1;
                    disp('gaussianMixBayesianLp: Using Martin"s edge-detector to initialize edges');
                else
                    flagEdgeDetectorFIX = 1;
                    disp('gaussianMixBayesianLp: and keep the fixed to that estimate!');
                end                
            otherwise
                fprintf('gaussianMixBayesianLp: Unknown option "%s" parsed\n', options(oCount));
        end
    end    
end
D = size(X, 1);
N = size(X, 2);
NhoodSize = 4;
NhoodDirections = 2; %WARNING:Unreliable constant

% Initialization
covar = zeros(D, D, K);
allDataCov = cov(X');
% Deterministic initialization -- K-means
[m w] = deterministicKmeans(X(:,1:N), K);
w = w'*ones(1, N);
w = w + 0.2*rand(K, N); disp('gaussianMixBayesianLp: Added peturbation noise on initialization.'); % Add some noise..% Add some noise..
w = (w ./ (ones(K, 1) * (sum(w, 1)+eps)))';
for i = 1:K
    % Make each std.deviation equal to 1/K of total std.deviation.
    covar(:,:,i) = allDataCov / K^2 + eps*eye(D);
end
%%% FIN - Deterministic initialization -- K-means %%%
beta = double(0.8 * ones(NhoodDirections, K));
u = double(ones(NhoodSize, K, N));
z = inv(N) * ones(N, K);
wDiff2 = zeros(NhoodSize, K, N);
wDiff2Beta = zeros(NhoodSize, K, N);
%%% Hyperparameter initialization
a0 = ones(NhoodSize, 1);
b0 = ones(NhoodSize, 1);
a_ksi = a0;
b_ksi = b0;

if flagEdgeDetector == 1
    for j = 1:K
        for n = 1:NhoodSize
            %% This will send [1,0] to [0, 1]
            u(n, j, :) = 1-edgemap(:);
        end
    end
end
% Variational inference iterations...
disp('VI iteration  Likelihood         Average Logl    LhoodIncrease   RandIndex   RandIndexContextual');
disp('------------------------------------------------------------------------------------------------');
likelihood = 0;
for iterations = 1:15
    prev = likelihood;
    likelihood = 0;    
    % E-step
    for i = 1:K
        z(:, i) = w(:, i) .* gaussianValue(X, m(:, i), covar(:,:,i)) + eps;
        likelihood = likelihood + z(:, i);
    end
    likelihood = sum(log(likelihood)); %logp(X|Pi), il nous faudra aussi +log(p(Pi))    
    z = z ./ (sum(z, 2) * ones(1, K));
    %% Special -- Compute the current segmentation score (Rand)
    if exist('groundTruth', 'var') && isempty(groundTruth) == 0
        [junk segmentation] = max(z'); [currentRI junk currentMCR] = randIndex(reshape(segmentation, imageSize), ...
            groundTruth);
        [junk segmentationCont] = max(w'); [currentRICont junk currentMCRCont] = randIndex(reshape(segmentationCont, imageSize), ...
            groundTruth);         
    else
        currentRI = nan; currentRICont = nan; currentMCR = nan; currentMCRCont = nan;
    end
    %% Special -- reinitialize w to the posterior
    if iterations == 1
        im=zeros(imageSize);
        for kk=1:K
            im(:,:) = reshape(z(:, kk), imageSize);
            xx=filter2((1./9.)*ones(3,3),im);
            z(:, kk)=xx(:);
        end
        w = z;
        %%%%
%         if exist('groundTruth', 'var') && isempty(groundTruth) == 0
%             disp('gaussianMix: Cheat! Initializing using ground truth!');
%             z = zeros(N, K);
%             for j = 1:K
%                 z(groundTruth == j, j) = 1;
%             end
%             w = z;
%         end
    end
    % M-step
    for i = 1:K
        if(sum(z(:,i)) > 0) %To catch empty clusters
            m(:, i) = (X * z(:, i)) / sum(z(:, i));
        end
        newCovar = inv(sum(z(:, i))) * ( ...
            (X - m(:, i)*ones(1,N))* ...
            sparse(1:N, 1:N, z(:, i), N, N) * ...
            (X - m(:, i)*ones(1,N))' );
        newCovar = newCovar + eps*eye(D);
        if rcond(newCovar) > 1e3 * eps
            covar(:, :, i) = newCovar;
        end
        clear newCovar;
    end
    %% Compute q*(U) %%
    % Some convenient statistics, regarding w. 
    %Neighbour square differences: wDiff2
    for n = 1:NhoodSize
        for j = 1:K
            [pos direction] = getNeighbourInfo(n);
            temp = translation(reshape(w(:, j), imageSize),  -pos);
            temp2 = translation(temp, +pos);    
            wDiff2(n, j, :) = (temp(:) - temp2(:)).^2;
            wDiff2Beta(n, j, :) = wDiff2(n, j, :) / (beta(direction, j) + eps);
            if (iterations > 1 || flagEdgeDetector == 0) && (flagEdgeDetectorFIX == 0)           
                u(n, j, :) = psi(a_ksi(n)) - psi(b_ksi(n));                        
                u(n, j, :) = squeeze(u(n, j, :)) + logGaussianValue((temp(:) - temp2(:))', 0, beta(direction, j));
                u(n, j, :) = sig(u(n, j, :));
            end
        end
    end
    %TODO: Likelihood (VLB) is far from complete
    llChangeRatio = (likelihood - prev) / abs(prev);
    fprintf('%3d           %3.2f         %2.5f        %2.5f%%        %2.4f       %2.4f        %2.4f       %2.4f\n', ...
        iterations, likelihood, likelihood/N, llChangeRatio*100, currentRI, currentRICont, ...
        currentMCR, currentMCRCont);
    % STEP 2: Maximize parameters: beta, v (freedom degrees) and pi
    %% Pi (weights) %%
    maxLevel = max([floor(log2(max(imageSize)/16)); 3]);    
    if flagGridScan == 0
        oldW = w;
        for j = 1:K
            aQuad = 0; bQuad = 0;
            cQuad = -0.5 * z(:, j); 
            for k = 1:NhoodSize
                [pos d] = getNeighbourInfo(k);
                temp = translation(reshape(oldW(:, j), imageSize),  -pos);
                aQuad = aQuad + inv(beta(d, j) + eps) * squeeze(u(k, j, :));
                bQuad = bQuad - inv(beta(d, j) + eps) * ...
                    (squeeze(u(k, j, :)) .* temp(:));
            end
            w(:, j) = solveQuad(aQuad, bQuad, cQuad);
        end
        for n = 1:N
            w(n, :) = BIDProjection(w(n, :));
        end
    else
        for mrLevel = maxLevel:-1:0
            w = MEXgridScan(mrLevel, imageSize, w, z, u, beta, K, flagGridScanAlsoZ);
        end        
    end
    %% Beta %%
    for j = 1:K
        betaComponent = zeros(1, NhoodSize);
        betaDenom = zeros(1, NhoodSize);
        for n = 1:NhoodSize
            betaComponent(n) = (sum(u(n, j, :) .* wDiff2(n, j, :)));
            betaDenom(n) = sum(u(n, j, :));
        end
        for d = 1:size(beta, 1)
            beta(d, j) = sum(betaComponent(getDirectionInfo(d))) / sum(betaDenom(getDirectionInfo(d)));
            if beta(d, j) < eps %Lowerbound for beta: 1e-10
                beta(d, j) = eps;
            end
        end
        clear betaComponent; clear betaDenom;
    end
%     for j = 1:K
%         betaComponent = zeros(1, NhoodSize);
%         for n = 1:NhoodSize
%             betaComponent(n) = (sum(u(n, j, :) .* wDiff2(n, j, :))) / ...
%                     getTotalNeighbours(getNeighbourInfo(n), imageSize);
%         end
%         for d = 1:size(beta, 1)
%             beta(d, j) = sum(betaComponent(getDirectionInfo(d)));
%         end
%         clear betaComponent;
%     end
    %% aksi, bksi %%
    for n = 1:NhoodSize
        a_ksi(n) = a0(n) + sum(u(n, :));
        b_ksi(n) = b0(n) + sum(1-u(n, :));
    end
    if(iterations == 5)
        %% Save the result for iterations=5.. this is for test purposes
        w5 = w';
        z5 = z';
    end            
end
w = w';
z = z';
return;

function res = getTotalNeighbours(offset, imageSize)
res = prod(imageSize - abs(offset));
return;

function [pos direction] = getNeighbourInfo(n)
switch n
    case {1}
        pos = [0 -1];
        direction = 1;
    case {2}
        pos = [-1 0];
        direction = 2;
    case {3}
        pos = [0 1];
        direction = 1;
    case {4}
        pos = [1 0];
        direction = 2;
    otherwise
        disp('Error: Unknown neighbour');
end
return;
function res = getDirectionInfo(d)
switch d
    case {1}
        res = [1 3];
    case {2}
        res = [2 4];
    otherwise
        disp('Error: Unknown direction');
end
return;

function res = solveQuad(a, b, c)
if a <= 0
    a = -a;
    b = -b;
    c = -c;
end
determining = b.^2 - 4 * a .* c;
res = (-b + sqrt(determining)) ./ (2*a+eps);
return;

function res = translation(rect, offset)
T = maketform('affine', [eye(2) [offset(1) offset(2)]'; 0 0 1]');
R = makeresampler('nearest', 'fill');
res = tformarray(rect, T, R, [1 2], [1 2], size(rect),[],[]);
return;

function res = sig(x)
res = 1./(1 + exp(-x));
return;