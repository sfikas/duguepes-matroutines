function [errorCode,message]=mex_aux(filename,option,option2);
% in case there is an #include file.cpp, put shared in the search
% path relative to location of mex file
% Timothee Cour, 29-Aug-2006 07:49:15


global data;
if isfield(data,'paths')
    isInclude=1;
    includeDir = data.paths.includeDir;
else
    isInclude=0;
    includeDir=pwd;
end

if nargout == 0
    isVerbose = 1;
else
    isVerbose = 0;
end

% % verbose:
% [errorCode,message]=mex_silent(isVerbose,'-v',['-I',includeDir],filename);%to display warnings 
% disp(message);
% return;

k=0;

k=k+1;args{k}=isVerbose;
k=k+1;args{k}='-O';
if ~isempty(option2)
    k=k+1;args{k}=option2;
end
if isInclude
    k=k+1;args{k}=['-I',includeDir];
end
k=k+1;args{k}=filename;
% if ~isempty(option)
%     k=k+1;args{k}=option;
% end
if ispc
    k=k+1;args{k}='-argcheck';
end

%% to export for matlab <7.1, use those 2 lines
if ispc
    k=k+1;args{k}=['-output'];
    k=k+1;args{k}=changeExt(filename,'.dll');
end

[errorCode,message]=mex_silent(args{:});
0;
%{
if isempty(option) && isempty(option2)
    [errorCode,message]=mex_silent(isVerbose,'-O',['-I',includeDir],filename);
elseif ~isempty(option) && isempty(option2)
    [errorCode,message]=mex_silent(isVerbose,'-O',['-I',includeDir],filename,option);
elseif isempty(option) && ~isempty(option2)
    [errorCode,message]=mex_silent(isVerbose,'-O',option2,['-I',includeDir],filename);
elseif ~isempty(option) && ~isempty(option2)
    [errorCode,message]=mex_silent(isVerbose,'-O',option2,['-I',includeDir],filename,option);
end
%}

% mex('-O',filename);
% mex(filename);


% if isempty(option) && isempty(option2)
%     mex('-O',['-I',includeDir],filename);
% elseif ~isempty(option) && isempty(option2)
%     mex('-O',['-I',includeDir],filename,option);
% elseif isempty(option) && ~isempty(option2)
%     mex('-O',option2,['-I',includeDir],filename);
% elseif ~isempty(option) && ~isempty(option2)
%     mex('-O',option2,['-I',includeDir],filename,option);
% end
