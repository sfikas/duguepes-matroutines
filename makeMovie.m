function res = makeMovie(SEG, timeDimension, colorMap);
% makeMovie
% Make a movie out of XxYxZ 'matrix' SEG.
% The 'Z' coordinate will be used as the
% time dimension in the movie.
%
% Example:
%       movie(makeMovie(A)); 
%   The following will make and play the movie for A using the y-dimension
%   as time, and a copper-tone color map:
%       movie(makeMovie(A, 2, copper));
%   The following will instead make the movie and save it externally
%   as an AVI file:
%       movie2avi(makeMovie(A), 'a_movie.avi');
% G.Sfikas 24 Mar 2007
%
if nargin < 3
    colorMap = YR_colors;
    if nargin < 2
        timeDimension = 3;
    end
end
if isempty(timeDimension)
    timeDimension = 3;
end
SEG = shiftdim(SEG, 3-timeDimension);
colormap(colorMap);
for slice = 1:size(SEG, 3)
imagesc(SEG(:,:,slice));
movieSEG(slice) = getframe();
end

res = movieSEG;
return;