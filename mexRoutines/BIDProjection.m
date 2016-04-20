function res = BIDProjection(x, w)
%
% BIDProjection(vector)
% Projects input vector onto constraint plane defined by
%  1. sum(vector) = 1
%  2. Each variate of 'vector' >= 0
%  vector must 1xK sized.
%
% If a 2nd argument 'weights' is defined, then the ConjugateProjection
% function is called. See also there for details.
% 
%  ***CET VERSION EST UN WRAPPER QUI APPELERA LA VERSION "MEX" DE CET FICHE.***
%  G.Sfikas 22 Mar 2008
%  based on code by K.Blekas
%
% Update 21 Apr 2008
if nargin == 1
%     res = MEX_BIDProjection(x, numel(x));
%     return;
    w = ones(size(x));
end
res = MEX_ConjugateProjection(x, w, numel(x));
return;