function [GCE]= GlobalConsistencyError(sampleLabels1,sampleLabels2)
% sampleLabels1 and sampleLabels2 must have 
% pixel values 1,2,3,...,NbClusters (must be segmentation maps)

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


% global consistancy error (from BSDS ICCV 01 paper) ... lower => better
N = sum(sum(n));
marginal_1 = sum(n,2);
marginal_2 = sum(n,1);
% the hackery is to protect against cases where some of the marginals are
% zero (should never happen, but seems to...)
E1 = 1 - sum( sum(n.*n,2) ./ (marginal_1 + (marginal_1 == 0)) ) / N;
E2 = 1 - sum( sum(n.*n,1) ./ (marginal_2 + (marginal_2 == 0)) ) / N;
GCE = min( E1, E2 );



