function [filelist filename] = getfilelist(directory,filter,nowarn)
% this function returns a matrix with each row is the address of one file,
% which is useful for passing to spm functions. 
%       FILELIST = getfilelist(DIRECTORY,FILTER)
% The files are selected within the DIRECTORY by FILTER which is a string 
% such as 'f*.img'.

if nargin<3, nowarn = 0; end;
if iscell(directory)
    m = length(directory);
    filelist = cell(m,1);
    filename = cell(m,1);
    for i=1:m
        [filelist{i},filename{i}]=getfilelist(directory{i},filter);
    end
    return;
end

currpath = pwd;
try
    cd(directory);
catch exception %#ok<NASGU>
    error('oops! cannot go inside: %s', directory);
end

list = dir(filter);
n = 0; 
for i=1:length(list)
    if(n<length(list(i).name))
        n = length(list(i).name);
    end
end
k = length(pwd) + 1;

%     filelist = char(length(list),k+n);
filelist = repmat(' ',length(list),k+n);
filename = repmat(' ',length(list),n);
for i=1:length(list)
    filelist(i,1:k+length(list(i).name)) = fullfile(pwd ,list(i).name);
    filename(i,1:length(list(i).name))   = list(i).name;
end
cd(currpath);
if(size(filelist,1)==0 && ~nowarn)
    warning('No file found!');
end

end
