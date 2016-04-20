function [res RI mcrMean MCR]= randIndex(sampleLabels1,sampleLabels2)
% [res RI mcrMean MCR] = randIndex(sampleLables1, sampleLabels2)
%
% Input:
% sampleLabels1 and sampleLabels2 must have 
% pixel values 1,2,3,...,NbClusters (must be segmentation maps)
% At most one of them, may include multiple segmentations,
% and can be a cell (like in the Berkely database).
%
% Output:
%       res         mean Rand index (over given ground truth maps)
%       RI          Rand indices
%       mcrMean     mean Misclassification ratio
%                   IMPORTANT: MCR computed is only approximate.
%                               The user is advised to doublecheck
%                               the output where possible;
%                               Otherwise, it may be better to use
%                               the Rand index only.
%       MCR         Misclassification ratios
%
% Example:
%   [randIndex junk mcrError] = randIndex(mySegmentation, myGroundTruth);
%
% Based on code by C.Nikou 10 Jul 2007
%
% TODO: Need some more work on handling 3-dimension images and more.
% G.Sfikas 7 Mar 2007.
% Update 15 Apr 2008, added code for misclassification ratio
% Update 23 Jan 2008, now can handle segmentations of different segments;
%                       correctly computing MCR.
% Handle multiple ground truths
%
% Update 13 Jan 2010, fixing bugs..(MCR)
nbGt = 1;
if iscell(sampleLabels2)
    tt = size(sampleLabels1);
    nbGt = numel(sampleLabels2);
    groundTruths = sampleLabels2;
else
    if numel(size(sampleLabels1)) > numel(size(sampleLabels2))
        tt = size(sampleLabels1);
        nbGt = tt(end);
        groundTruths = sampleLabels1; sampleLabels1 = sampleLabels2;
    elseif numel(size(sampleLabels1)) < numel(size(sampleLabels2))
        tt = size(sampleLabels2);
        nbGt = tt(end);
        groundTruths = sampleLabels2;
    end
end

RI = zeros(1, nbGt); MCR = zeros(1, nbGt);
for jj = 1:nbGt
    if nbGt > 1
        if numel(tt) == 2
            if iscell(groundTruths) == 1
                sampleLabels2 = groundTruths{jj};
            else
                sampleLabels2 = groundTruths(:, jj);
            end
        elseif numel(tt) == 3
            if iscell(groundTruths) == 1
                sampleLabels2 = groundTruths{jj};
            else        
                sampleLabels2 = groundTruths(:, :, jj);
            end
        else
            error('randIndex: Cannot handle ground truth size.');
        end
    end
    [imWidth,imHeight]=size(sampleLabels1);
    if (numel(sampleLabels1) == numel(sampleLabels2))
        sampleLabels2 = reshape(sampleLabels2, size(sampleLabels1));
    end
    [imWidth2,imHeight2]=size(sampleLabels2);
    N=imWidth*imHeight;
    if (imWidth~=imWidth2)||(imHeight~=imHeight2)
        sampleLabels2 = imresize(sampleLabels2, size(sampleLabels1), 'nearest');
    end;

    % make the group indices start at 1
    if min(min(sampleLabels1)) < 1
        sampleLabels1 = sampleLabels1 - min(min(sampleLabels1)) + 1;
    end
    if min(min(sampleLabels2)) < 1
        sampleLabels2 = sampleLabels2 - min(min(sampleLabels2)) + 1;
    end

    segmentcount1=max(max(sampleLabels1));
    segmentcount2=max(max(sampleLabels2));

    % compute the count matrix
    %  from this we can quickly compute rand index, GCE, VOI, ect...
    n=zeros(segmentcount1,segmentcount2);

    for i=1:imWidth
        for j=1:imHeight
            u=sampleLabels1(i,j);
            v=sampleLabels2(i,j);
            n(u,v)=n(u,v)+1;
        end;
    end;


    % the rand index, in [0,1] ... higher => better
    %  fast computation is based on equation (2.2) of Rand's paper.
    N = sum(sum(n));
    n_u=sum(n,2);
    n_v=sum(n,1);
    N_choose_2=N*(N-1)/2;
    RI(jj) = 1 - ( sum(n_u .* n_u)/2 + sum(n_v .* n_v)/2 - sum(sum(n.*n)) )/N_choose_2;
    % The misclassification ratio, in [0,1] ... lower = better!
    % This is a fast (greedy) 'approximation'.
    % For the real MCR we need to compute all permutations of the
    %       'n' matrix.
    correctGuys1 = 0; n1 = n;
    for i = 1:size(n, 1) %Search by every row
        [Y I] = max(n1(i, :));
        correctGuys1 = correctGuys1 + Y;
        n1(:, I(1)) = zeros(size(n1, 1), 1);
    end
    correctGuys2 = 0; n2 = n;
    for i = 1:size(n, 2) %Search by every column
        [Y I] = max(n2(:, i));
        correctGuys2 = correctGuys2 + Y;
        n2(I(1), :) = zeros(1, size(n2, 2));
    end    
    MCR(jj) = 1 - max(correctGuys1, correctGuys2)/N;
end
res = mean(RI);
mcrMean = mean(MCR);
return;