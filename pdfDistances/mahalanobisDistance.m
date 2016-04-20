function res = mahalanobisDistance(descQuery, descSet, blobWorldBinCorelations);
%
%   Calculate distance of a blob to an image (set of blobs).
%   Using distance described in BlobWorld paper (Quadratic-form)
%   Note that the result is x'Ax , _not_ x'inv(A)x (ie 'mahalanobis' is
%       a misnomer for that matter)
%
%   G.Sfikas 14/6/2006
%   Update 28/6/2006 (no compound query)
%   Update/fix 24/10/2006
%   Update comments 7/2/2007
if nargin < 3
    %if 'bWBC' argument _is_ supplied, things run a bit faster.
    blobWorldBinCorelations = getBlobWorldBinCorelations();
end
dif = (ones(size(descSet,1), 1) * descQuery) - descSet;
d1 = sum(dif * blobWorldBinCorelations .* dif, 2);
res = min(d1);
return;