function y = scaleSelection(inimage);
% Select 'proper' smoothing scale (see Blobworld)
% using polarity measure.
% 'inimage' contains L* values.
% Output contains standard deviation values,
%  corresponding to proper gaussian convolution kernels
%  for each pixel.
%
% y = scaleSelection(inimage)
%
% G.Sfikas 13/4/2006
% 
% fix: 02/05/2006
%
threshold = 0.02; %As specified in Blobworld paper
for k = 1:8
    scale = (k-1)/2;
    A = convolution2D(inimage, scale);
    [tpl, junk, junk2] = computePolarity(inimage, scale);
    polarity(:,:,k) = convolution2D(tpl, 2*scale);
end
polarityGrad = diff(polarity, 1, 3);
% stopMap and stopIndex contain info about the correct stopping scale.
stopMap = abs(polarityGrad) <= threshold;
stopMap(:,:,end+1) = 1;
for m = 1:size(stopMap,3)
    % The second addent is a trick to disregard zeros on the 'min'
    % operation later.
    stopMap2(:,:,m) = stopMap(:,:,m)*m + (stopMap(:,:,m)==0)*(size(stopMap,3)+1);
end
[temp, stopIndex] = min(stopMap2, [], 3);
stopIndex = (stopIndex-1)/2;

y = stopIndex;
return;

