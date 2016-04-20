function pr = probabilisticRandIndex(Im,ImBerkeley)
% pr = probabilisticRandIndex(Im, ImBerkeley)
% Computes the PR index between segmentation Im and
% ground truths in ImBerkeley. Higher index is better.
% Segment indices must be in the range [1..+Inf].
%
% Im            Y x X segmentation (K-level image).
% ImBerkeley    NbGT x Y x X. That is NbGT ground truths, where their
%                   segment numbers may vary.
% pr            pr index, a score in [0, 1].
%
% Adaptation of the function courtesy of C.Nikou. (See student-t directory)
%
% See also
%       randIndex
%
% G. Sfikas 16 Nov 2007

disp('probabilisticRandIndex: The present code for this function is obsolete(too slow).');
disp('probabilisticRandIndex: You are advised to use "randIndex" instead.');
if numel(size(ImBerkeley)) < 3
    disp('You must provide more than one ground truths.');
    disp('For index for a single ground truth, see "randIndex".');
end
NbGT = size(ImBerkeley, 1);
[R C]=size(Im);

%Store the number of clusters of each ground truth segmentation map
GTClusters=zeros(NbGT,1);
for k=1:NbGT
    x=ImBerkeley(k,:,:);
    GTClusters(k)=max(x(:));
end

pr=0;
for k=1:NbGT
    s=0;
    for i=1:R
        for j=1:C


            for m=i+1:R
                for n=1:C

                    if Im(i,j)==Im(m,n)
                        c=1;
                    else
                        c=0;
                    end

                    if ImBerkeley(k,i,j)==ImBerkeley(k,m,n)
                        I1=1;
                    else
                        I1=0;
                    end

                    I2=1-I1;

                    s = s + c*I1+(1-c)*I2;
                end
            end

            m=i;
            for n=1:C-1
                if Im(i,j)==Im(m,n)
                    c=1;
                else
                    c=0;
                end

                if ImBerkeley(k,i,j)==ImBerkeley(k,m,n)
                    I1=1;
                else
                    I1=0;
                end

                I2=1-I1;

                s = s + c*I1+(1-c)*I2;

            end


        end
    end
    pr = pr + s;
end

pr = pr / (NbGT*(R*C)*(R*C-1)/2.0);