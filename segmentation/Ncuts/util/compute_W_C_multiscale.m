function [W,C]=compute_W_C_multiscale(image);
% compute multiscale image affinity matrix W and multiscale constraint
% matrix C from input image
% Timothee Cour, 29-Aug-2006 07:49:15


[p,q,r] = size(image);
dataW = computeParametersW(image);
layers = computeParametersLayers(p,q);

[C,C12]=computeMultiscaleConstraints(layers);

% compute each layers(i).location as subsamples of the finest layer
layers=computeLocationFromConstraints(C12,layers);

W=computeMultiscaleW(image,layers,dataW);
