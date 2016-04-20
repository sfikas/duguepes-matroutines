function y = convertJxN(A)
% Convert a MxNxJ matrix to JxMN. The conversion is done columnwise.
% Can handle matrices of higher dimensionality, eg 
% can convert MxNxPxJ to JxMNP.
% Example:
%       convertJxN(imread('lenna.jpg'))
%   will give a 3xMN matrix, w/ each column the rgb intensities of a
%   single pixel.
%
% G.Sfikas 29/5/2006
%
% Update 16 Mar 2007: Handles matrices of any dimensionality too, now.
% (Eg 4-D, 5-D etc).
% Fix 16 Mar 2007: y = shiftdim(t2(:,:)', +1) should be t2(:,:) (no
% transpose). _DOUBTFUL_
t2 = shiftdim(A, max(size(size(A))) - 1);
y = shiftdim(t2(:,:),+1)';
return;