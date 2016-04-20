function res = buildBlobDescriptor(labColor, contrast, anisotropy, means, vars);
% Build a blob descriptor as described in blobworld paper.
%
% Example:
%   buildBlobDescriptor(labColor, contrast, anisotropy).
%
% -------------------------------------------------------------------------
% labColor:     3xN matrix, w/ each row containing L*a*b* values for 
%               a single pixel.
% contrast:     1xN or Nx1 vector containing contrast values.
% anisotropy:   1xN or Nx1 vector containing anisotropy values.
% Returns:      A 502x1 vector (the descriptor).
% -------------------------------------------------------------------------
%
% See also
%
% G.Sfikas 13/11/06
%
%   UPDATE HISTORY:
%   ---------------
%   13-11-06 Symmazema
%   22-11-06 Allages gia normalization

if isempty(labColor)
    res = zeros(1,502);
    return;
end
res = hist(findProperBin(labColor), 0:499);
res = res / (sum(res));
% Normalization
contrastMean = means(6);
contrastVar = vars(6);
caMean = means(9);
caVar = vars(9);
temp1 = (contrast - contrastMean) / contrastVar;
res(501) = mean(temp1);
temp = contrast .* anisotropy;
temp = (temp - caMean) / caVar;
res(502) = mean(temp);
return;