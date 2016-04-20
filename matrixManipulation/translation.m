function res = translation(rect, offset, filltype)
% res = translation(rect, offset, filltype)
%
% Translate (move) the matrix 'rect' _plus_ the 'offset'.
%
% rect                          Input matrix.
% offset                        Vector of directions to move 'rect'.
% filltype      {'fill'}        Zero-padding.
%               'replicate'     Padding using border values.
%
% Examples
%       translation([1 2 3; 4 5 6; 7 8 9], [1 0]) will produce
%               [0 0 0;
%                1 2 3;
%                4 5 6]
%   while
%       translation([1 2 3; 4 5 6; 7 8 9], [1 0],'replicate') will pad like
%               [1 2 3;
%                1 2 3;
%                4 5 6]
%
% G.Sfikas 15 Feb 2008
if nargin < 3
    filltype = 'fill';
end
T = maketform('affine', [eye(2) [offset(1) offset(2)]'; 0 0 1]');
R = makeresampler('nearest', filltype);
res = tformarray(rect, T, R, [1 2], [1 2], size(rect),[],[]);
return;
