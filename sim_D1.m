function sim_D1(i)
% effects of number of trials

mode = '';
if nargin==0, mode = 'sumfig'; end
if nargin==1, mode = 'run'; end
switch mode
    case 'run'
        run(i);
    case 'sumfig'
        sumfig;
    otherwise
        error('Unknown mode: %s',mode);        
end
end

function run(ii)
nsim        = 20;
if nargin<1, ii=nsim; end
nsim        = 20;

simcat      = 'sim_D';
Nbar        = [10 30];

symname     = 'sim_rw';

T1          = 50;
simstr1     = sprintf('%s%sT%d',symname,mat2str(Nbar),T1);

T2          = 200;
simstr2     = sprintf('%s%sT%d',symname,mat2str(Nbar),T2);

models = {'model_rw1','model_rw2'};
tasks  = {'task_rw1','task_rw2'};
pnames = { {'\alpha','\beta'},{'\alpha^+','\alpha^-','\beta'} };

modelnames = {'Regular RL','Dual-\alpha RL'}; 
mnames = modelnames;

tau1 = 4*[1 1]';
tau2 = 4*[1 1 1]';
tau  = {tau1,tau2};

tx1  = [.1 1]';
tx2  = [.8 .4 3]';
logit = @(x)log(x./(1-x));

mu1  = [logit(tx1(1)); log(tx1(2))];
mu2  = [logit(tx2(1:2)); log(tx2(3))];

mu   = {mu1,mu2};

normx1 = {@(x)1./(1+exp(-x)), @exp};
normx2 = {@(x)1./(1+exp(-x)),@(x)1./(1+exp(-x)),@exp};

normx  = {normx1,normx2};

sim_sim(ii,simcat,simstr1,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx,T1);
if ii(end)~=nsim, return; end
sim_wrap(nsim,simcat,simstr1,modelnames,mnames,pnames,normx);

sim_sim(ii,simcat,simstr2,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx,T2);
if ii(end)~=nsim, return; end
sim_wrap(nsim,simcat,simstr2,modelnames,mnames,pnames,normx);
end

function sumfig
fsum    = fullfile(getdefaults('sumdir'),'sim_D','sum.mat');

if ~exist(fsum,'file')
    simcats = {'sim_D','sim_C','sim_D'};
    simstrs = {'fit_sim_rw[10 30]T50','fit_sim_rw[10 30]','fit_sim_rw[10 30]T200'};

    fsimfits = fullfile(getdefaults('pipedir'),simcats,simstrs);

    [bms, mdx, sdx, methods] = sim_stat_D1(fsimfits);
    xgroups = {'T=50','T=100','T=200'};
    save(fsum,'xgroups','bms','mdx','sdx','methods');
else
    sums    = load(fsum);
    xgroups = sums.xgroups;
    bms     = sums.bms;
    mdx     = sums.mdx;
    sdx     = sums.sdx;
    methods = sums.methods;    
end

bms(bms<1) = 1;

% close all;
fig(xgroups,bms,mdx,sdx,methods);
set(gcf,'name',mfilename);
end

function fig(xgroups,bms,mdx,sdx,methods)
[fs,fst,fsl,fsy,fsalpha,fsxt,fsA,fn,fnt,xsA,ysA,fpos0,siz0,colmap,alf,cmaphbi,bw,colmapsim]=fig_plot_properties; %#ok<*ASGLU>

ysA = 1.2;

fpos = [0.4 .25];
siz = fpos./fpos0(3:4).*siz0;

% bw  = .27/.68*fpos0(3);
bwi = .75*bw;

nr = 1;
nc = 2;

figure;
set(gcf,'units','centimeters');
fsiz = get(gcf,'position');
fsiz(3:4) = siz;
fsiz(1:2) = fpos0(1:2).*fsiz(1:2);
set(gcf,'position',fsiz);
set(gcf,'units','normalized');

%---------
subplot(nr,nc,1);    
errorbarKxN(bms,bms*0,xgroups,'',colmap,0,bwi);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
title(sprintf('Model selection\nperformance'),'fontsize',fst,'fontname',fnt);
ylim([0 110]);
ytick = 0:20:100;
set(gca,'ytick',ytick);
ylabel('Correct selection %','fontsize',fsy);
text(xsA,ysA,'A','fontsize',fsA,'Unit','normalized','fontname','Calibri');
ax = ancestor(gca, 'axes');
xaxes = get(ax,'XAxis');
set(xaxes,'fontsize',fsxt);
%--------
subplot(nr,nc,2);    
errorbarKxN(mdx,sdx,xgroups,'',colmap,0,bwi);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
title(sprintf('Parameter estimation\nperformance'),'fontsize',fst,'fontname',fnt);
% ylim([0 1.1]);
ytick = 0:.2:1;
set(gca,'ytick',ytick);
ylabel('Error','fontsize',fsy);
text(xsA,ysA,'B','fontsize',fsA,'Unit','normalized','fontname','Calibri');
ax = ancestor(gca, 'axes');
xaxes = get(ax,'XAxis');
set(xaxes,'fontsize',fsxt);

lg = legend(methods,'location','northeast','fontsize',fsl,'Orientation','horizontal');
set(lg,'linewidth',1);

end