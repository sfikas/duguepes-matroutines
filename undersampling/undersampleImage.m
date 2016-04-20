function res = undersampleImage(inimage, desiredDim1, desiredDim2)
%
% undersampleImage(inimage, desiredDim1, desiredDim2);
%
A = inimage;
dim1 = size(A,1);
dim2 = size(A,2);
if (dim1 < dim2 && desiredDim1 > desiredDim2) || ...
        (dim1 > dim2 && desiredDim1 < desiredDim2)
    tempD = desiredDim2;
    desiredDim2 = desiredDim1;
    desiredDim1 = tempD;
end
maskX = ones(dim1,1)*undersamplingMask(dim2, desiredDim2);
maskY = undersamplingMask(dim1, desiredDim1)'*ones(1,dim2);
finalMask = maskX & maskY;
finalMask = finalMask(:);

channel1 = A(:,:,1);
channel1 = channel1(finalMask);
channel1 = reshape(channel1, desiredDim1, desiredDim2);
channel2 = A(:,:,2);
channel2 = channel2(finalMask);
channel2 = reshape(channel2, desiredDim1, desiredDim2);
channel3 = A(:,:,3);
channel3 = channel3(finalMask);
channel3 = reshape(channel3, desiredDim1, desiredDim2);

res = zeros(desiredDim1, desiredDim2, 3);
res(:,:,1) = channel1;
res(:,:,2) = channel2;
res(:,:,3) = channel3;
res = uint8(res);
return;



function res = undersamplingMask(origsize, newsize)
A = 1:origsize;
res = diff(floor((1:origsize) / (origsize/newsize)));
if(size(find(res), 2) < newsize)
    res(origsize) = 1;
else
    res(origsize) = 0;
end

return;