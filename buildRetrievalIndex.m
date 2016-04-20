function res = buildRetrievalIndex(gmmdir, extension)
% Creates an index of all (typically, model) filenames.
% Uses information from gmmdir directory.
%
% Examples:
%       buildRetrievalIndex('gmm', '.gmm')
%   or
%       buildRetrievalIndex()
%   will give a char matrix w/ each row containing the filename of
%   a single gaussian mixture model.
%       buildRetrievalIndex(fullfile(matlabroot,'work'),'.m')
%   will give a char matrix w/ each row containing the filename of
%   a MATLAB 'm' file on the program's 'work' directory.
%
% See also doRetrieval.
%
%   UPDATE HISTORY
%   --------------
%   08-05-06
%   15-06-06
%   14-11-06 symmazema / update
%
% G.Sfikas 8-5-2006
%
if nargin < 2
    extension = '.gmm';
    if nargin < 1
        gmmdir = 'gmm';
    end
end
savedir = pwd;
% t = clock;
index = char;       % Afto kanei ton index keno pinaka.
dir_struct = dir(gmmdir);
cd(gmmdir);
dirsize = size(dir_struct, 1);
for i = 1:dirsize
    if (dir_struct(i).isdir == 1) & (strcmp(dir_struct(i).name,'.') ~= 1) & (strcmp(dir_struct(i).name,'..') ~= 1)
        % Kolla ta arxeia tou epomenou katalogou sth lista index.
        t1 = dig(dir_struct(i).name, extension);
        index2 = concat(index, t1);
        clear index;
        index = index2;
    end
end
res = index;
cd(savedir);
% etime(clock,t)
return

function resdig = dig(directory, extension)
j = 1;
resdig = char;
cd(directory);
dir_struct = dir(pwd);
dirsize = size(dir_struct, 1);
for i = 1:dirsize
    [path, name, ext, ver] = fileparts(dir_struct(i).name);
    if (dir_struct(i).isdir ~= 1) & (strcmp(ext,extension) == 1)
        %
        % To padding to vazw gia na exoun ola ta strings
        % idio mhkos- etsi tha mpoyn ston pinaka 'resdig'
        % xwris provlima.
        clear tmpresdig;
        tmpresdig(1,:) = fullfile(pwd, dir_struct(i).name);
        fnlen = size( fullfile(pwd, dir_struct(i).name), 2);
        for padding = 1:(150 - fnlen)
            tmpresdig(1,fnlen+padding) = ' ';
        end
        resdig(j,:) = tmpresdig(1,:);
        j = j + 1;
    end
end
cd('..');
return

function resconcat = concat(list1, list2)
size1 = size(list1, 1);
size2 = size(list2, 1);
resconcat = list1;
resconcat((size1 + 1):(size1 + size2),1:size(list2,2)) = list2(1:size2,1:size(list2,2));
return
