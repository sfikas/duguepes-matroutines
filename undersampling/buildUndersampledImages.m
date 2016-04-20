function buildUndersampledImages(desiredDim1, desiredDim2);
% Builds undersampled versions of images in the 'pix'
% directory, and saves them on 'und'
%
% eg buildUndersampledImages(128,192);
%
% All images will be scaled so that:
%
% Dim1 -> desiredDim1
% Dim2 -> desiredDim2
%
% Where Dim1 = max(Dim1, Dim2);
%       Dim2 = min(Dim1, Dim2);
%
% Update 17/1/2006: erased 'cd to root' line.
%
pixdir = 'pix';
savedir = pwd;
tic
index = char;       % Afto kanei ton index keno pinaka.
dir_struct = dir(pixdir);
cd(pixdir);
dirsize = size(dir_struct, 1);
for i = 1:dirsize
    if (dir_struct(i).isdir == 1) & (strcmp(dir_struct(i).name,'.') ~= 1) & (strcmp(dir_struct(i).name,'..') ~= 1)        
        dig(dir_struct(i).name, desiredDim1, desiredDim2);
        disp(sprintf('Completed %s', dir_struct(i).name));
    end
end
cd(savedir);
toc
return

function resdig = dig(directory, desiredDim1, desiredDim2)
cd(directory);
dir_struct = dir(pwd);
dirsize = size(dir_struct, 1);
for i = 1:dirsize
    [path, name, ext, ver] = fileparts(dir_struct(i).name);
    if (dir_struct(i).isdir ~= 1) & ((strcmp(ext,'.jpg') == 1) | (strcmp(ext,'.jpeg') == 1) ...
            | strcmp(ext,'.JPG') == 1 | strcmp(ext,'.JPEG') == 1)
        inputFile = fullfile(pwd, dir_struct(i).name);
        undersampledData = undersampleImage(imread(inputFile), desiredDim1, desiredDim2);        
        imwrite(undersampledData, image2underimage(inputFile), 'jpeg');
        disp(sprintf('Done with image %s.', inputFile));
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