function [meanScores res] = diceIndex(sampleLabels1,sampleLabels2)
% sampleLabels1 and sampleLabels2 must have 
% pixel values 1,2,3,...,NbClusters (must be segmentation maps)
%
% See also
%       randIndex
% 
% G.Sfikas 28 Feb 2008

[imWidth,imHeight]=size(sampleLabels1);
[imWidth2,imHeight2]=size(sampleLabels2);
N=imWidth*imHeight;
if (imWidth~=imWidth2)||(imHeight~=imHeight2)
    disp( 'Input sizes: ' );
    disp( size(sampleLabels1) );
    disp( size(sampleLabels2) );
    error('Input sizes do not match in compare_segmentations.m');
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Until here, the code was the same as 'randIndex' !
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(n, 1)
    V1 = sum(n(i, :));
    V2 = sum(n(:, i));
    V12 = n(i, i);
    res(i) = 2*V12/(V1 + V2);
end
res = res(:);

meanScores = mean(res);
return;