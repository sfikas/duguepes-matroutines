function y = smoothUsingVariantScale(inimage, scaleMatrix)
% Smooths image 'inimage' using a different smoothing scale
%  for each pixel. The proper standard deviation values
%  must be supplied in matrix 'scaleMatrix'.
% 'inimage' are L* or a or b values.
% 'scaleMatrix' contains values up to 3.5.
%
% y = smoothUsingVariantScale(inimage, scaleMatrix)
%
% G.Sfikas 02/05/2006
%
scaleMatrix2 = scaleMatrix*2 + 1;
finalImage = inimage.*(scaleMatrix2==1);
% scale '0' is no smoothing.
for u = 2:8
    scale = (u-1)/2;
    finalImage = finalImage + convolution2D(inimage, scale).*(scaleMatrix2==u);
end
y = finalImage;
return;