function [VI]= VariationOfInformation(sampleLabels1,sampleLabels2)
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



% variation of information a "distance", in (0,vi_max] ... lower => better
N = sum(sum(n));
joint = n / N; % the joint pmf of the two labels
marginal_2 = sum(joint,1);  % row vector
marginal_1 = sum(joint,2);  % column vector
H1 = - sum( marginal_1 .* log2(marginal_1 + (marginal_1 == 0) ) ); % entropy of the first label
H2 = - sum( marginal_2 .* log2(marginal_2 + (marginal_2 == 0) ) ); % entropy of the second label
MI = sum(sum( joint .* log2_quotient( joint, marginal_1*marginal_2 )  )); % mutual information
VI = H1 + H2 - 2 * MI; 


% log2_quotient 
%   helper function for computing the mutual information
%   returns a matrix whose ij entry is
%     log2( a_ij / b_ij )   if a_ij, b_ij > 0
%     0                     if a_ij is 0
%     log2( a_ij + 1 )      if a_ij > 0 but b_ij = 0 (this behavior should
%                                         not be encountered in practice!
function lq = log2_quotient( A, B )
lq = log2( (A + ((A==0).*B) + (B==0)) ./ (B + (B==0)) );
