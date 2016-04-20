function image = imread2(filename,maxSize,isGray);
% Timothee Cour, 29-Aug-2006 07:49:15

image = imread(filename);
if nargin >=3 && isGray
    if size(image,3)>1
        image = rgb2gray(image);
    end
end


if nargin > 1 && ~isempty(maxSize)
    [p,q,r] = size(image);
    if max(p,q)>maxSize
        factor = maxSize/max(p,q);
        image = imresize(image,round([p,q]*factor),'bicubic');
    end
end

image = rescaleImage(double(image));
