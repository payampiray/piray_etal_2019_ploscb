function [fs,fst,fsl,fsy,fsalpha,fsxt,fsA,fn,fnt,xsA,ysA,fpos0,siz,colmap,alf,cmaphbi,bw,colmapsim]=fig_plot_properties
fsbase  = 13;
fd = 5;

fs      = fsbase;
fst     = fsbase+fd;
fsl     = fsbase+fd;
fsy     = fsbase+fd;
fsalpha = fsbase+fd+2;
fsxt    = fsbase+fd-2; % xticklabel
fsA     = fsbase+fd;

% fs      = 8;
% fst     = 8;
% fsl     = 8;
% fsy     = 8;
% fsalpha = 8;
% fsxt    = 8; % xticklabel
% fsA     = 8;


% fn  = 'MyriadPro-Regular';
fn  = 'Calibri';
fnt = 'calibri';

xsA = -.15;
ysA = 1.1;

% fpos0 = [0    0.0800    .68    0.63];
fpos0 = [.2    0.2    .68    0.63];
siz = [46 24]; %cm

colmap = [.2 0 0; .8 .4 .1; 1 .2 .2] ;
alf = .6;
cmaphbi = colmap(3,:);

bw  = .27/.68*fpos0(3);

colmapsim = [.2 .6 1;.05 .15 .3;0 .4 .65; .7 .85 .95];
end