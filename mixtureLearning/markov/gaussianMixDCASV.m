function [m, covar, w, z, u] = gaussianMixDCASV(X, K, imageSize, groundTruth)
% [m, covar, w, z, u] = gaussianMixDCASV(X, maxSegments,
%                           imageSize 
%       | 2 |
%     --/   \--
%     1   x   3     <-- Neighbourhood indexing
%     --\   /--
%       | 4 |
% G.Sfikas 26 Nov 2007
% Updates 01 Feb 2007

D = size(X, 1);
N = size(X, 2);
NhoodSize = 4;
NhoodDirections = 2; %WARNING:Unreliable constant

% Initialization
[m w junk covar] = deterministicKmeans(X(:,1:N), K);
w = w'*ones(1, N);
w = (w ./ (ones(K, 1) * (sum(w, 1)+eps)))';
beta = double(0.8 * ones(NhoodDirections, K));
u = double(ones(NhoodSize, K, N));
%%% Hyperparameter initialization
z = inv(N) * ones(N, K);
% newWeight = inv(N) * ones(N, K);
wDiff2 = zeros(NhoodSize, K, N);


% disp('gaussian Mix: on labels.'); 
% w = inv(K) * ones(size(w));
% EM iterations...
disp('EM iteration  Likelihood         Average Logl    LhoodIncrease   RandIndex   RandIndexContextual');
disp('------------------------------------------------------------------------------------------------');
likelihood = 0;
for iterations = 1:15
    if mod(iterations, 229) == 0
        fprintf('.');
    end
    prev = likelihood;
    likelihood = 0;
    % E-step
    for i = 1:K
        z(:, i) = w(:, i) .* gaussianValue(X, m(:, i), covar(:,:,i));
        likelihood = likelihood + z(:, i);
    end
    likelihood = sum(log(likelihood)); %logp(X|Pi), il nous faudra aussi +log(p(Pi))
    z = z ./ (sum(z, 2) * ones(1, K));
    %% Special -- Compute the current segmentation score (Rand)
    if exist('groundTruth', 'var') && isempty(groundTruth) == 0
        [junk segmentation] = max(z'); currentRI = randIndex(segmentation, groundTruth(:)');
        [junk segmentationCont] = max(w'); currentRICont = randIndex(segmentationCont, groundTruth(:)'); 
    else
        currentRI = nan; currentRICont = nan;
    end
    %% Special -- reinitialize w to the posterior
    if iterations == 1
%         im=zeros(imageSize);
%         for kk=1:K
%             im(:,:) = reshape(z(:, kk), imageSize);
%             xx=filter2((1./9.)*ones(3,3),im);
%             z(:, kk)=xx(:);
%         end
%         w = z;
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
    % Some convenient statistics, regarding w. 
    %Neighbour square differences: wDiff2
    for j = 1:K
        for n = 1:NhoodSize
            [pos direction] = getNeighbourInfo(n);
            temp = translation(reshape(w(:, j), imageSize),  -pos);
            temp2 = translation(temp, +pos);
            wDiff2(n, j, :) = (temp(:) - temp2(:)).^2;
            likelihood = likelihood + ...
                sum(logGaussianValue((temp(:) - temp2(:))', 0, beta(direction, j)));
        end
    end
    llChangeRatio = (likelihood - prev) / abs(prev);
    fprintf('%3d           %3.2f         %2.5f        %2.5f%%        %2.4f       %2.4f\n', ...
        iterations, likelihood, likelihood/N, llChangeRatio*100, currentRI, currentRICont);
    % STEP 2: Maximize parameters: beta, v (freedom degrees) and pi
    %% Pi (weights) %%
    for ii = 1:2
    for j = 1:K
        aQuad = 0; bQuad = 0;
        cQuad = -0.5 * z(:, j); 
        for k = 1:NhoodSize
            [pos d] = getNeighbourInfo(k);
            temp = translation(reshape(w(:, j), imageSize),  -pos);
            aQuad = aQuad + inv(beta(d, j) + eps) * squeeze(u(k, j, :));
            bQuad = bQuad - inv(beta(d, j) + eps) * ...
                (squeeze(u(k, j, :)) .* temp(:));
        end
        newWeight(:, j) = solveQuad(aQuad, bQuad, cQuad);
    end
    for n = 1:N
        w(n, :) = BIDProjection(newWeight(n, :));
    end
    end
%     for n1 = 1:imageSize(1)
%     for n2 = 1:imageSize(2)
%         n = (n2-1)*imageSize(1) + n1;
%         for j = 1:K
%             aQuad = 0; bQuad = 0;
%             cQuad = -0.5 * z(n, j);
%             for k = 1:NhoodSize
%                 [pos d] = getNeighbourInfo(k);
%                 neighPos = [n1+pos(1) n2+pos(2)];
%                 if neighPos(1) < 1 || neighPos(1) > imageSize(1) || ...
%                         neighPos(2) < 1 || neighPos(2) > imageSize(2)
%                     continue;
%                 end
%                 neighPos = (neighPos(2)-1)*imageSize(1) + neighPos(1);
%                 aQuad = aQuad + inv(beta(d, j) + eps) * squeeze(u(k, j, n));
%                 bQuad = bQuad - inv(beta(d, j) + eps) * ...
%                     (squeeze(u(k, j, n)) .* w(neighPos, j));
%             end
%             newWeight(j) = solveQuad(aQuad, bQuad, cQuad);
%         end
%         w(n, :) = BIDProjection(newWeight);
%     end
%     end
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

function res = solveQuad(a, b, c)
if a < 0
    a = -a;
    b = -b;
    c = -c;
end
determining = b.^2 - 4 * a .* c;
res = (-b + sqrt(determining)) ./ (2*a);
return;