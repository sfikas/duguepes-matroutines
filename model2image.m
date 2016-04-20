function imagename = model2image(modelname, K)
% Convert a modelname into the corresponding imagename, etc
% D:\MATLAB7\work\main\gmm\nonbear\1254.gmm
%  converts to
% D:\MATLAB7\work\main\pix\nonbear\1254.jpg
%
% H KANEI to akrivws antistrofo. Analogws.
%
% Update 20/6/2006
% 
% Update 23 Nov 2007 (K argument)
% G.Sfikas
imagepath = 'pix';
modelpath = 'gmm';
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
if nargin > 1
    part3 = sprintf('%s%d', part3, K);
end
part1 = modelname(1:(lastind2-1));
part4 = modelname(ind+1:size(modelname, 2));
if (strcmp(ext(1:4), '.seg') == 1) || (strcmp(ext(1:4), '.gmm') == 1)
    imagename = fullfile(part1, imagepath, part3, strcat(name,'.jpg'));
else
    imagename = fullfile(part1, modelpath, part3, strcat(name,'.gmm'));
end    
return