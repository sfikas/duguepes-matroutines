function [m, covar, w, z, u, v, beta, likelihood, w5, z5] = gaussianMixBayesianContinuousLp(X, K, imageSize, groundTruth, options, edgemap)
% [m, covar, w, z, u] = gaussianMixBayesianContinuousLp(X, maxSegments,
%                           imageSize
% 
% Gaussian mixture model fit to data X. A markov random field is assumed on the
% spatially-varying weights, with a _Student-t_ field clique function.
% 4-neighbourhood is considered on the MRF. 
%       | 2 |
%     --/   \--
%     1   x   3     <-- Neighbourhood indexing
%     --\   /--
%       | 4 |
% G.Sfikas 20 Nov 2007
%
% Update history
%       1. 13 May 2008
%       2. 25 Oct 2008 (warning message for adding noise on init)
if (exist('edgemap', 'var') == 0)
    disp('gaussianMixBayesianContinuousLp: Will work without an initial edgemap');
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
                disp('gaussianMixBayesianContinuousLp: Using grid-scan optimization on Pi');                
            case {'z'}
                flagGridScanAlsoZ = 1;
                disp('gaussianMixBayesianContinuousLp: Using grid-scan also on Z');                
            case {'m'}
                flagMultiRes = 1;
                disp('gaussianMixBayesianContinuousLp: Using multiresolution optimization (DEPRECATED)');
            case {'e'}
                if flagEdgeDetector == 0
                    flagEdgeDetector = 1;
                    disp('gaussianMixBayesianContinuousLp: Using Martin"s edge-detector to initialize edges');
                else
                    flagEdgeDetectorFIX = 1;
                    disp('gaussianMixBayesianContinuousLp: and keep the fixed to that estimate!');
                end
            otherwise
                fprintf('gaussianMixBayesianContinuousLp: Unknown option "%s" parsed\n', options(oCount));
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
w = w + 0.2*rand(K, N); disp('gaussianMixBayesianContinuousLp: Added peturbation noise on initialization.'); % Add some noise..
w = (w ./ (ones(K, 1) * (sum(w, 1)+eps)))';
for i = 1:K
    % Make each std.deviation equal to 1/K of total std.deviation.
    covar(:,:,i) = allDataCov / K^2 + eps*eye(D);
end
%%% FIN - Deterministic initialization -- K-means %%%
beta = double(0.8 * ones(NhoodDirections, K));
u = double(ones(NhoodSize, K, N));
logU = zeros(NhoodSize, K, N);
v = double(ones(NhoodDirections, K));
%%% Hyperparameter initialization
z = inv(N) * ones(N, K);
wDiff2 = zeros(NhoodSize, K, N);
wDiff2Beta = zeros(NhoodSize, K, N);


if flagEdgeDetector == 1
    for j = 1:K
        for n = 1:NhoodSize
            %% This will send [1,0] to (0, +inf)
            u(n, j, :) = -log(edgemap(:) + eps) + 2*eps;
            logU(n, j, :) = log( u(n,j,:));
        end
    end
end
% EM iterations...
disp('EM iteration  Likelihood   Average Logl   LhoodIncrease   RandIndex');
disp('---------------------------------------------------------------------');
likelihood = 0;
for iterations = 1:15 % Normalement cela est fixe a quinze
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
        [junk segmentation] = max(z'); [currentRI junk currentMCR] = randIndex(reshape(segmentation, imageSize), groundTruth);
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
    for j = 1:K
        for n = 1:NhoodSize
            [pos direction] = getNeighbourInfo(n);
            temp = translation(reshape(w(:, j), imageSize),  -pos);
%             temp2 = translation(reshape(w(:, j), imageSize), +pos);
            temp2 = translation(temp, +pos);
            wDiff2(n, j, :) = (temp(:) - temp2(:)).^2;
            wDiff2Beta(n, j, :) = wDiff2(n, j, :) / (beta(direction, j) + eps);
            likelihood = likelihood + ...
                sum(logStudentValue((temp(:) - temp2(:))', 0, beta(direction, j) + eps, v(direction, j)));
        end
    end
    llChangeRatio = (likelihood - prev) / abs(prev);
    fprintf('%3d           %3.2f         %2.5f        %2.5f%%        %2.4f       %2.4f        %2.4f       %2.4f\n', ...
        iterations, likelihood, likelihood/N, llChangeRatio*100, currentRI, currentRICont, ...
        currentMCR, currentMCRCont);
    if (iterations > 1 || flagEdgeDetector == 0) && (flagEdgeDetectorFIX == 0)
        for j = 1:K
            for n = 1:NhoodSize
                u(n, j, :) = (v(j) + 1) ./ (v(j) + wDiff2Beta(n, j, :));
                logU(n, j, :) = ...
                    psi((v(j) + 1)/2) - log((v(j) + wDiff2Beta(n, j, :))/2);
            end
        end
    end
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
        for n = 1:NhoodSize
            betaComponent(n) = (sum(u(n, j, :) .* wDiff2(n, j, :))) / ...
                    getTotalNeighbours(getNeighbourInfo(n), imageSize);
        end
        for d = 1:size(beta, 1)
            beta(d, j) = sum(betaComponent(getDirectionInfo(d)));
        end
        clear betaComponent;
    end
    %% V (Degrees of freedom) %%
    for j = 1:K
        vComponent = zeros(1, NhoodSize);
        for n = 1:NhoodSize
            vComponent(n) = (sum(logU(n, j, :) - u(n, j, :)) / ...
                getTotalNeighbours(getNeighbourInfo(n), imageSize));
        end
        for d = 1:size(v, 1)
            vConstant = sum(vComponent(getDirectionInfo(d))) + 1;
            v(d, j) = fzero(@(vTemp) vFunction(vTemp, vConstant),[+eps +inv(eps)]);
        end
        clear vComponent;
    end
    % STEP 3: Evaluate the lower bound.
    %% TODO .. plus tard ..
    if(iterations == 5)
        %% Save the result for iterations=5.. this is for test purposes
        w5 = w';
        z5 = z';
    end            
end
w = w';
z = z';
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

function res = getTotalNeighbours(offset, imageSize)
res = prod(imageSize - abs(offset));
return;

function res = vFunction(v, vConstant)
res = log(0.5*v) - psi(0.5*v) + vConstant;
return;

function res = solveQuad(a, b, c)
if a < 0
    a = -a;
    b = -b;
    c = -c;
end
determining = b.^2 - 4 * a .* c;
res = (-b + sqrt(determining)) ./ (2*a);
return;

function res = translation(rect, offset)
T = maketform('affine', [eye(2) [offset(1) offset(2)]'; 0 0 1]');
R = makeresampler('nearest', 'fill');
res = tformarray(rect, T, R, [1 2], [1 2], size(rect),[],[]);
return;