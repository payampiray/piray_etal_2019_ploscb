function varargout = getdefaults(action,varargin)

st = dbstack('-completenames');
maindir = (fileparts(st(1).file));
pipedir = fullfile(maindir,'pipe');
tempdir = pipedir;
sumdir  = fullfile(maindir,'sum');

switch action
    case 'welcome'
        clc; close all;
        st = dbstack;
        [~,mfname] = fileparts(st(2).file);
        fprintf('%-40s%30s\n',mfname,datestr(now));
        fprintf('%-70s\n\n',repmat('=',1,70));
    case 'pipedir'
        varargout{:} = pipedir;
    case 'tempdir'
        varargout{:} = tempdir;
    case 'sumdir'
        varargout{:} = sumdir;
    case 'addpath'
        addpath(fullfile(maindir,'tools'));
    otherwise
        error('%s does not exist',action);
end
    

end
