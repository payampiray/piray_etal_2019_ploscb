function sim_G1(i)
% Figure 8
% big model space (4 models)

mode = '';
if nargin==0, mode = 'sumfig'; end
if nargin==1, mode = 'run'; end
getdefaults('addpath');

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
for j=1:4
Nbar = [10 10 10 10];
Nbar(j) = 30;
run_Nbar(Nbar,ii,nsim);
end

end

function run_Nbar(Nbar,ii,nsim)

simcat      = 'sim_G';

symname     = 'big';
simstr      = sprintf('%s%s',symname,mat2str(Nbar));

models     = {'model_rw1','model_rw2','model_kf1','model_ac1'};
tasks      = {'task_rw1' ,'task_rw2','task_kf1','task_ac1'};
modelnames = {'Regarlar RL','Dual-\alpha RL','Kalman filter','Actor-Critic'};
mnames     = modelnames;

pnames = { {'\alpha','\beta'},...
           {'\alpha^+','\alpha^-','\beta'} ,...
           {'\omega','\beta'} ,...
           {'\alpha^{\it c}','\alpha^{\it a}','\beta'}};

tau1 = 4*[1 1]';
tau2 = 4*[1 1 1]';
tau3 = 4*[1 1]';
tau4 = 4*[1 1 1]';
tau  = {tau1,tau2,tau3,tau4};

tx1  = [.1 1]';
tx2  = [.8 .4 3]';
tx3  = [2 2]';
tx4  = [.1 .6 3]';

logit  = @(x)log(x./(1-x));
mu1    = [logit(tx1(1)); log(tx1(2));];
mu2    = [logit(tx2(1:2)); log(tx2(3));];
mu3    = log(tx3(1:2));
mu4    = [logit(tx4(1:2)); log(tx4(3));];

mu     = {mu1,mu2,mu3,mu4};

normx1 = {@(x)1./(1+exp(-x)),@exp};
normx2 = {@(x)1./(1+exp(-x)),@(x)1./(1+exp(-x)),@exp};
normx3 = {@exp,@exp};
normx4 = {@(x)1./(1+exp(-x)),@(x)1./(1+exp(-x)),@exp};

normx  = {normx1,normx2,normx3,normx4};

sim_sim(ii,simcat,simstr,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx);
if ii(end)~=nsim, return; end
sim_wrap(nsim,simcat,simstr,modelnames,mnames,pnames,normx);
end

function sumfig
simcat = 'sim_G';
fsum    = fullfile(getdefaults('sumdir'),simcat,'sum_G1.mat');
if ~exist(fsum,'file')
    fitsimstrs = {'fit_big[30 10 10 10]','fit_big[10 30 10 10]','fit_big[10 10 30 10]','fit_big[10 10 10 30]'}; 
    fsimfits  = fullfile(getdefaults('pipedir'),simcat,fitsimstrs);
    [mpxp,spxp,mnbar,snbar,mbms,mdmx,sdmx, ~,~,methods] = sim_stat_G1(fsimfits);
    mnames = {'RL','DA','KF','AC'};
    simnames = {'Scenario 1','Scenario 2','Scenario 3','Scenario 4'};
    save(fsum,'mpxp','spxp','mnbar','snbar','mbms','mdmx','sdmx','methods','mnames','simnames');
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
    simnames  = sums.simnames;
end

% close all;
fig(mpxp,spxp,mnbar,snbar,mbms,mdmx,sdmx,methods,mnames,simnames);
set(gcf,'name',mfilename);
end

function fig(mpxp,spxp,mnbar,snbar,mbms,mdmx,sdmx,methods,mnames,simnames)

[fs,fst,fsl,fsy,fsalpha,fsxt,fsA,fn,fnt,xsA,ysA,fpos0,siz0,colmap,alf,cmaphbi,bw,colmapsim]=fig_plot_properties; %#ok<*ASGLU>

fpos = [0.6 fpos0(4)];
siz = fpos./fpos0(3:4).*siz0;

bwi = .75*bw;
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
errorbarKxN(mpxp,spxp,simnames,'',colmapsim,0,bwi);
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
errorbarKxN(mnbar,snbar,simnames,'',colmapsim,0,bwi);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
title('Model comparison by HBI','fontsize',fst,'fontname',fnt);
ylabel('Model frequency','fontsize',fsy);
text(xsA,ysA,'B','fontsize',fsA,'Unit','normalized','fontname','Calibri');
ax = ancestor(gca, 'axes');
xaxes = get(ax,'XAxis');
set(xaxes,'fontsize',fsxt);

legend(mnames,'location','north','fontsize',fsl,'Orientation','horizontal');

yl = get(gca,'ylim');
ytick = get(gca,'ytick'); dy = ytick(2)-ytick(1); ytick(end) = [];
ylim([0 yl(2)+dy*.5]);
set(gca,'ytick',ytick);
%--------
subplot(nr,nc,3);    
errorbarKxN(mbms,0*mbms,simnames,'',colmap,0,bwi);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
title('Model selection performance','fontsize',fst,'fontname',fnt);
ylabel('Correct selection %','fontsize',fsy);
ax = ancestor(gca, 'axes');
xaxes = get(ax,'XAxis');
set(xaxes,'fontsize',fsxt);

text(xsA,ysA,'C','fontsize',fsA,'Unit','normalized','fontname','Calibri');

ylim([0 105]); ytick = 0:20:110;
set(gca,'ytick',ytick);
%--------
subplot(nr,nc,4 );    
errorbarKxN(mdmx,sdmx,simnames,'',colmap,0,bwi);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
title('Parameter estimation performance','fontsize',fst,'fontname',fnt);
ylabel('Error','fontsize',fsy);
ax = ancestor(gca, 'axes');
xaxes = get(ax,'XAxis');
set(xaxes,'fontsize',fsxt);

legend(methods,'location','north','fontsize',fsl,'Orientation','horizontal');
text(xsA,ysA,'D','fontsize',fsA,'Unit','normalized','fontname','Calibri');

myl = 1.01*max(max(mdmx+sdmx));
ylim([0 myl]);

end