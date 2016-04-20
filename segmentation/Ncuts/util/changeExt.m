function path2=changeExt(path1,ext2);
% Timothee Cour, 29-Aug-2006 07:49:15

if ~strcmp(ext2(1),'.')
    error('ext2 must start with .');
end
[path,name,ext]=fileparts(path1);
path2=fullfile(path,[name,ext2]);

