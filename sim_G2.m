function sim_G2(i,cat)
% Figure 9
% a more complicated task: 2-step task

mode = '';
if nargin==0, mode = 'sumfig'; end
if nargin==2, mode = 'run'; end
getdefaults('addpath');

switch mode
    case 'run'
        run(i,cat);
    case 'sumfig'
        sumfig;
    otherwise
        error('Unknown mode: %s',mode);        
end
end

function run(ii,cat)
nsim        = 20;
if nargin<1, ii=nsim; end
nsim        = 20;

simcat      = 'sim_G';
Nbar        = [30 10 10];

symname     = '2step';
simstr      = sprintf('%s%s',symname,mat2str(Nbar));

models     = {'model_hybrid','model_mb','model_mf'};
tasks      = {'task_hybrid','task_mb','task_mf'};

fxnames    = {'\alpha_1', '\alpha_2', '\lambda', '\beta_1', '\it w', '\it p', '\beta_2'};
pnames     = { fxnames , fxnames([2 4 6 7]), fxnames([1 2 3 4 6 7]) };

modelnames = {'hyrbid','model-based','model-free'}; 
mnames     = {'Hyrbid','MB','MF'};

tau1 = 4*ones(7,1);
tau2 = 4*ones(4,1);
tau3 = 4*ones(6,1);
tau  = {tau1,tau2,tau3};

tx  = [0.5 0.4 0.6 5.2 .4 .1 3.7];
logit = @(x)log(x./(1-x));
mu1 = [logit(tx(1:3)) log(tx(4)) logit(tx(5)) tx(6) log(tx(7))]';
mu2 = mu1([2 4 6 7]);
mu3 = mu1([1 2 3 4 6 7]);

mu  = {mu1,mu2,mu3};

normx1 = {@(x)1./(1+exp(-x)), @(x)1./(1+exp(-x)), @(x)1./(1+exp(-x)), @exp, @(x)1./(1+exp(-x)), @(x)x, @exp};
normx2 = normx1([2 4 6 7]);
normx3 = normx1([1 2 3 4 6 7]);

normx  = {normx1, normx2, normx3};

frandwalks = {'2step_rndwalk1.mat','2step_rndwalk2.mat'};
frandwalks = fullfile(getdefaults('pipedir'),simcat,frandwalks);

simcat0 = 'sim_2step';
simstr0 = 'daw3o[30 10 10]';

sim_sim_G2(ii,simcat,simstr,frandwalks,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx,simcat0,simstr0,cat);
if ii(end)~=nsim && cat==8
sim_wrap(nsim,simcat,simstr,modelnames,mnames,pnames,normx);
end

end

function sumfig
simcat = 'sim_G';
fsum    = fullfile(getdefaults('sumdir'),simcat,'sum_G2.mat');
if ~exist(fsum,'file')
    fitsimstrs = {'fit_2step[30 10 10]'};
    fsimfits  = fullfile(getdefaults('pipedir'),simcat,fitsimstrs);
    [mpxp,spxp,mnbar,snbar,mbms,~,~, mdmx,sdmx,methods] = sim_stat_G1(fsimfits);
    mnames     = {'Hybrid','MB','MF'};
    mnbar =[mnbar'; .6 .2 .2];
    snbar =[snbar'; 0 0 0];
    Nbarstrs = {'HBI','True'};
    save(fsum,'mpxp','spxp','mnbar','snbar','mbms','mdmx','sdmx','methods','mnames','Nbarstrs');    
else
    sums    = load(fsum);
    mpxp    = sums.mpxp;
    spxp    = sums.spxp;
    mnbar   = sums.mnbar;
    snbar   = sums.snbar;
    mbms    = sums.mbms;
    mdmx    = sums.mdmx;
    sdmx    = sums.sdmx;
    methods = sums.methods;
    mnames  = sums.mnames;
    Nbarstrs  = sums.Nbarstrs;    
end

% close all;
fig(mpxp,spxp,mnbar,snbar,Nbarstrs,mbms, mdmx,sdmx,methods,mnames);
set(gcf,'name',mfilename);
end

function fig(mpxp,spxp,mnbar,snbar,Nbarstrs,mbms,mdx,sdx,methods,mnames)
[fs,fst,fsl,fsy,fsalpha,fsxt,fsA,fn,fnt,xsA,ysA,fpos0,siz0,colmap,alf,cmaphbi,bw,colmapsim]=fig_plot_properties; %#ok<*ASGLU>
cmaphbi2 = [cmaphbi; colmapsim(2,:)];

xsA = -.2;
ysr = .1;

fpos = [0.4 fpos0(4)];
siz = fpos./fpos0(3:4).*siz0;

bwi = 1.2*bw;

%---------

nr = 2;
nc = 2;

figure;
set(gcf,'units','centimeters');
fsiz = get(gcf,'position');
fsiz(3:4) = siz;
fsiz(1:2) = fpos0(1:2).*fsiz(1:2);
set(gcf,'position',fsiz);
set(gcf,'units','normalized');

subplot(nr,nc,1);    
errorbarKxN(mpxp',spxp',mnames,'',cmaphbi,0,bwi);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
title('Model selection by HBI','fontsize',fst,'fontname',fnt);
ylim([0 1.1]);
ytick = 0:.2:1;
set(gca,'ytick',ytick);
ylabel('PXP','fontsize',fsy);
text(xsA,ysA,'A','fontsize',fsA,'Unit','normalized','fontname','Calibri');

ax = ancestor(gca, 'axes');
xaxes = get(ax,'XAxis');
set(xaxes,'fontsize',fsxt);
%--------
subplot(nr,nc,2 );    
errorbarKxN(mnbar,snbar,mnames,'',cmaphbi2,0,bwi);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
title('Model comparison by HBI','fontsize',fst,'fontname',fnt);
ylabel('Model frequency','fontsize',fsy);
text(xsA,ysA,'B','fontsize',fsA,'Unit','normalized','fontname','Calibri');

yl = get(gca,'ylim');
ylim([yl(1) yl(2)-.01]);

ax = ancestor(gca, 'axes');
xaxes = get(ax,'XAxis');
set(xaxes,'fontsize',fsxt);

lg = legend(Nbarstrs,'location','northeast','fontsize',fsl,'Orientation','vertical'); hold on;
set(lg,'linewidth',1);
%--------
subplot(nr,nc,3);    
errorbarNxK(mbms',0*mbms',methods,'',colmap,0,bwi);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
title(sprintf('Model selection\nperformance'),'fontsize',fst,'fontname',fnt);
ylabel('Correct selection %','fontsize',fsy);

ylim([0 105]); ytick = 0:20:110;
set(gca,'ytick',ytick);

ax = ancestor(gca, 'axes');
xaxes = get(ax,'XAxis');
set(xaxes,'fontsize',fsxt);

text(xsA,ysA+ysr,'C','fontsize',fsA,'Unit','normalized','fontname','Calibri');
%--------
mx = mdx{1}(:,5);
ex = sdx{1}(:,5);

subplot(nr,nc,4 );    
errorbarNxK(mx',ex',methods,'',colmap,0,bwi);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
title(sprintf('Parameter estimation \nof the weight parameter'),'fontsize',fst,'fontname',fnt);
ylabel('Error','fontsize',fsy);

ax = ancestor(gca, 'axes');
xaxes = get(ax,'XAxis');
set(xaxes,'fontsize',fsxt);

text(xsA,ysA+ysr,'D','fontsize',fsA,'Unit','normalized','fontname','Calibri');

end
