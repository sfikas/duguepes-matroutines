function y = convolution2D(A, scale)
% Convolution w/ Gaussian kernel of standard deviation 'scale'
%
% y = convolution2D(scale)
%
% G.Sfikas 13/5/2006
%
%
kernel = gaussKernel1D(scale);
y = conv2(conv2(A, kernel, 'same'), kernel', 'same');
return;

function y = gaussKernel1D(sdeviation)
% Returns a maxSizex1 Gaussian kernel, normalized to sum 1,
% with standard deviation 'sdeviation'.
%
% y = gaussKernel1D(sdeviation)
%
% G.Sfikas 11/4/2006
%
maxSize = 13;
m = ceil(maxSize/2);
if sdeviation <= 0
    kernel = zeros(maxSize, 1);
    kernel(m) = 1;
    y = kernel;
    return;
end
for x = 1:maxSize
    kernel(x) = exp(-0.5*sdeviation^-2*(x-m)^2);
end
kernel = kernel';
kernel = kernel/sum(kernel);
y = kernel;

return;