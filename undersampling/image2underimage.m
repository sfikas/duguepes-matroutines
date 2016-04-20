function underimagename = image2underimage(modelname)
% Convert an imagename into the corresponding undersampled imagename, etc
% D:\MATLAB7\work\main\pix\nonbear\1254.jpg
%  converts to
% D:\MATLAB7\work\main\und\nonbear\1254.jpg
%
% Update 29/6/2006
%
imagepath = 'pix';
underpath = 'und';
[path,name,ext,ver] = fileparts(modelname);
ind = 1;
lastind = 1;
lastind2 = 1;
for i = 1:size(modelname, 2);
    if modelname(i) == '\' | modelname(i) == '/'
        lastind2 = lastind;
        lastind = ind;
        ind = i;
    end
end
part3 = modelname(lastind+1:ind-1);
part1 = modelname(1:(lastind2-1));
dirname = modelname((lastind2+1):(lastind-1));
part4 = modelname(ind+1:size(modelname, 2));
if (strcmp(dirname, 'pix') == 1)
    underimagename = fullfile(part1, underpath, part3, strcat(name,'.jpg'));
else
    underimagename = fullfile(part1, imagepath, part3, strcat(name,'.jpg'));
end    
return