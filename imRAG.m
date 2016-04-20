function varargout = imRAG(img, varargin)
%IMRAG compute region adjacency graph of labeled image
%
%   Usage :
%   ADJ = imRAG(IMG);
%   compute region adjacencies graph of labeled image IMG. The result is
%   a N*2 array, containing 2 indices for each couple of neighbor regions.
%   ADJ as the format [LBL1 LBL2], LBL1 and LBL2 being vertical arrays the
%   same size.
%
%   LBL1 is given in ascending order, LBL2 is given in ascending order for
%   each LBL1. Ex :
%   [1 2]
%   [1 3]
%   [1 4]
%   [2 3]
%   [2 5]
%   [3 4]
%
%   [NODES, ADJ] = imRAG(IMG);
%   Return two arrays: the first one is a [N*2] array containing centroids
%   of the N labeled region, and ADJ is the adjacency previously described.
%   For 3D images, the nodes array is [N*3].
%   
%   Requires the image processing toolbox for computing centroids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on the work of David Legland (INRIA)
% 
% The difference is that this version does not need to find
% a black pixel in the between-superpixels area to count two
% pixels as neighbours. I do that automatically.
%
%

edgeMap = edge(img, 'sobel', 0);
img(edgeMap) = 0;

% count number of regions
N = double(max(img(:)));

% size of image
dim = size(img);

% initialize array
edges = [];

if length(dim)==2
	% compute matrix of absolute differences in the first direction
	diff1 = abs(diff(double(img), 1, 1));
	
	% find non zero values (region changes)
	[i1 i2] = find(diff1);
	
	% delete values close to border
	i2 = i2(i1<dim(1)-1);
	i1 = i1(i1<dim(1)-1);
	
	% get values of consecutive changes
	val1 = diff1(sub2ind(size(diff1), i1, i2));
	val2 = diff1(sub2ind(size(diff1), i1+1, i2));
	
	% find changes separated with 2 pixels
	ind = find(val2 & val1~=val2);
	edges = unique([val1(ind) val2(ind)], 'rows');	
	
	
	% compute matrix of absolute differences in the second direction
	diff2 = abs(diff(double(img), 1, 2));
	
	% find non zero values (region changes)
	[i1 i2] = find(diff2);
	
	% delete values close to border
	i1 = i1(i2<dim(2)-1);
	i2 = i2(i2<dim(2)-1);
	
	% get values of consecutive changes
	val1 = diff2(sub2ind(size(diff2), i1, i2));
	val2 = diff2(sub2ind(size(diff2), i1, i2+1));
	
	% find changes separated with 2 pixels
	ind = find(val2 & val1~=val2);
	edges = [edges ; unique([val1(ind) val2(ind)], 'rows')];
	
elseif length(dim)==3
    
     % compute matrix of absolute differences in the first direction
	diff1 = abs(diff(double(img), 1, 1));
	
	% find non zero values (region changes)
	[i1 i2 i3] = ind2sub(size(diff1), find(diff1));
	
	% delete values close to border
	i2 = i2(i1<dim(1)-1);
	i3 = i3(i1<dim(1)-1);
	i1 = i1(i1<dim(1)-1);
	
	% get values of consecutive changes
	val1 = diff1(sub2ind(size(diff1), i1, i2, i3));
	val2 = diff1(sub2ind(size(diff1), i1+1, i2, i3));
	
	% find changes separated with 2 pixels
	ind = find(val2 & val1~=val2);
	edges = unique([val1(ind) val2(ind)], 'rows');	
	
	
	% compute matrix of absolute differences in the second direction
	diff2 = abs(diff(double(img), 1, 2));
	
	% find non zero values (region changes)
	[i1 i2 i3] = ind2sub(size(diff2), find(diff2));
	
	% delete values close to border
	i1 = i1(i2<dim(2)-1);
	i3 = i3(i2<dim(2)-1);
	i2 = i2(i2<dim(2)-1);
	
	% get values of consecutive changes
	val1 = diff2(sub2ind(size(diff2), i1, i2, i3));
	val2 = diff2(sub2ind(size(diff2), i1, i2+1, i3));
	
	% find changes separated with 2 pixels
	ind = find(val2 & val1~=val2);
	edges = [edges ; unique([val1(ind) val2(ind)], 'rows')];	
	
	
	% compute matrix of absolute differences in the third direction
	diff3 = abs(diff(double(img), 1, 3));
	
	% find non zero values (region changes)
	[i1 i2 i3] = ind2sub(size(diff3), find(diff3));
	
	% delete values close to border
	i1 = i1(i3<dim(3)-1);
	i2 = i2(i3<dim(3)-1);
	i3 = i3(i3<dim(3)-1);
	
	% get values of consecutive changes
	val1 = diff3(sub2ind(size(diff3), i1, i2, i3));
	val2 = diff3(sub2ind(size(diff3), i1, i2, i3+1));
	
	% find changes separated with 2 pixels
	ind = find(val2 & val1~=val2);
	edges = [edges ; unique([val1(ind) val2(ind)], 'rows')];
end

% format output to have increasing order of n1,  n1<n2, and
% increasing order of n2 for n1=constant.
edges = sortrows(sort(edges, 2));

% remove eventual double edges
edges = unique(edges, 'rows');


if nargout == 1
    varargout{1} = edges;
end

if nargout == 2
    % Also compute region centroids
    stats = regionprops(img, 'centroid');
    tab = [stats.Centroid];
    if length(size(img))==2
        points = [tab(1:2:2*N-1)' tab(2:2:2*N)'];
    elseif length(size(img))==3
        points = [tab(1:3:3*N-2)' tab(2:3:3*N-1)' tab(3:3:3*N)'];
    end
    
    
    varargout{1} = points;
    varargout{2} = edges;
end
