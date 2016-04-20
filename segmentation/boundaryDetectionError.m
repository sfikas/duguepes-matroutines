function [BDE] = BoundaryDetectionError(imageLabels1,imageLabels2);

[imageX, imageY] = size(imageLabels1);
if imageX~=size(imageLabels2,1) | imageY~=size(imageLabels2,2)
    error('The sizes of the two comparing images must be the same.');
end

if isempty(find(imageLabels1~=imageLabels1(1)))
    % imageLabels1 only has one group
    boundary1 = zeros(size(imageLabels1));
    boundary1(1,:) = 1;
    boundary1(:,1) = 1;
    boundary1(end,:) = 1;
    boundary1(:,end) = 1;
else
    % Generate boundary maps
    [cx,cy] = gradient(imageLabels1);
    [boundaryPixelX{1},boundaryPixelY{1}] = find((abs(cx)+abs(cy))~=0);
    
    boundary1 = abs(cx) + abs(cy) > 0;
end

if isempty(find(imageLabels2~=imageLabels2(1)))
    % imageLabels2 only has one group
    boundary2 = zeros(size(imageLabels2));
    boundary2(1,:) = 1;
    boundary2(:,1) = 1;
    boundary2(end,:) = 1;
    boundary2(:,end) = 1;    
else    
    % Generate boundary maps
    [cx,cy] = gradient(imageLabels2);
    [boundaryPixelX{2},boundaryPixelY{2}] = find((abs(cx)+abs(cy))~=0);
    
    boundary2 = abs(cx) + abs(cy) > 0;
end

% boundary1 and boundary2 are now binary boundary masks. compute their
% distance transforms:
D1 = bwdist(boundary1);
D2 = bwdist(boundary2);

% compute the distance of the pixels in boundary1 to the nearest pixel in
% boundary2:
dist_12 = sum(sum(boundary1 .* D2 ));
dist_21 = sum(sum(boundary2 .* D1 ));

avgError_12 = dist_12 / sum(sum(boundary1));
avgError_21 = dist_21 / sum(sum(boundary2));

BDE = (avgError_12 + avgError_21) / 2;
