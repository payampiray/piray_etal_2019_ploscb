function sim_E1(i)
% Figure 6
% effects of number of subjects

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
nsim        = 500;
if nargin<1, ii=nsim; end

run_N([3 9],ii,nsim);
run_N([4 12],ii,nsim);
run_N([5 15],ii,nsim);
run_N([6 18],ii,nsim);

run_N([9 3],ii,nsim);
run_N([12 4],ii,nsim);
run_N([15 5],ii,nsim);
run_N([18 6],ii,nsim);
end

function run_N(Nbar,ii,nsim)

simcat      = 'sim_E';
symname     = 'sim_rw';
simstr      = sprintf('%s%s',symname,mat2str(Nbar));


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


sim_sim(ii,simcat,simstr,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx);
if ii(end)~=nsim, return; end
sim_wrap(nsim,simcat,simstr,modelnames,mnames,pnames,normx);
end

function sumfig
simcat   = 'sim_E';
fitsimstrs(1,:) = {'fit_sim_rw[3 9]','fit_sim_rw[4 12]','fit_sim_rw[5 15]',...
                   'fit_sim_rw[6 18]'};
fitsimstrs(2,:) = {'fit_sim_rw[9 3]','fit_sim_rw[12 4]','fit_sim_rw[15 5]',...
                   'fit_sim_rw[18 6]'};
xgroups = 12:4:24;
ratio   = .75;

fsum  = fullfile(getdefaults('sumdir'),simcat,sprintf('N%s.mat',mat2str(xgroups)));
if ~exist(fsum,'file')
    fsimfits  = fullfile(getdefaults('pipedir'),simcat,fitsimstrs);   
    
    [mpxp, elpxp, ehpxp, mnbar, elnbar, ehnbar, ms95, ms50, mdx, eldx, ehdx, methods_rfx, X, Y, auc, methods]= sim_stat_E1(fsimfits);
    save(fsum,'mpxp','elpxp','ehpxp','mnbar', 'elnbar', 'ehnbar', ...
              'ms95','ms50', 'mdx', 'eldx', 'ehdx','methods_rfx','X','Y','auc','methods');
else
    sums = load(fsum);
    mpxp        = sums.mpxp;
    elpxp       = sums.elpxp;
    ehpxp       = sums.ehpxp;
    mnbar        = sums.mnbar;
    elnbar       = sums.elnbar;
    ehnbar       = sums.ehnbar;
    ms95          = sums.ms95;
    ms50          = sums.ms50;
    methods_rfx = sums.methods_rfx;

    mdx        = sums.mdx;
    eldx       = sums.eldx;
    ehdx       = sums.ehdx;
    auc         = sums.auc;
    methods     = sums.methods;
end

ms95(ms95<1) = 1;
ms50(ms50<1) = 1;

% close all;
fig_plot_E(ratio,xgroups,mpxp, elpxp, ehpxp,mnbar, elnbar, ehnbar, ms95, ms50, methods_rfx, mdx, eldx, ehdx, auc, methods);
set(gcf,'name',mfilename);
end

function fig_plot_E(ratio,xgroups,mpxp, elpxp, ehpxp,mnbar, elnbar, ehnbar, ms95, ms50, methods_rfx, mdx, eldx, ehdx, auc, methods)
[fs,fst,fsl,fsy,fsalpha,fsxt,fsA,fn,fnt,xsA,ysA,fpos0,siz0,colmap,alf,cmaphbi,bw,colmapsim]=fig_plot_properties; %#ok<*ASGLU>

fpos = [.4 .85];
siz = fpos./fpos0(3:4).*siz0;

% ysA = 1.2;
ysr = .1;
ysA = ysA+ysr;

% bwi = .6*bw;
bwe = .9*bw;


figure;
set(gcf,'units','centimeters');
fsiz = get(gcf,'position');
fsiz(3:4) = siz;
fsiz(1:2) = fpos0(1:2).*fsiz(1:2);
set(gcf,'position',fsiz);
set(gcf,'units','normalized');
%---------

nr = 3;
nc = 2;
% N  = length(xgroups);

subplot(nr,nc,1);
errorbarKxN(mpxp,[elpxp; ehpxp],xgroups,'',colmap([1 3],:),0,bw);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
title('Exceedance probability','fontsize',fst,'fontname',fnt);
ylim([0 1.1]);
ytick = 0:.2:1;
set(gca,'ytick',ytick);
ylabel('PXP','fontsize',fsy);
text(xsA,ysA,'A','fontsize',fsA,'Unit','normalized','fontname','Calibri');

%------------
subplot(nr,nc,2);
errorbarKxN(mnbar,[elnbar; ehnbar],xgroups,'',colmap([1 3],:),0,bw);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
title('Model comparison','fontsize',fst,'fontname',fnt);
ylabel('Model frequency','fontsize',fsy);
ylim([0 1.1]);
ytick = 0:.2:1;
set(gca,'ytick',ytick);

hold on;
xtick = get(gca,'xlim');
plot(xtick,ratio*ones(size(xtick)),'color','k');

text(xsA,ysA,'B','fontsize',fsA,'Unit','normalized','fontname','Calibri');
%------------

subplot(nr,nc,3);
errorbarKxN(ms50,ms50*0,xgroups,'',colmap([1 3],:),0,bw);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
title('Model selection at PXP>0.5','fontsize',fst,'fontname',fnt);
ylabel(sprintf('Correct selection %%'),'fontsize',fsy);
ylim([0 115]);
ytick = 0:20:100;
set(gca,'ytick',ytick); 

text(xsA,ysA,'C','fontsize',fsA,'Unit','normalized','fontname','Calibri');

%------------

subplot(nr,nc,4);
errorbarKxN(ms95,ms95*0,xgroups,'',colmap([1 3],:),0,bw);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
title('Model selection at PXP>0.95','fontsize',fst,'fontname',fnt);
ylabel(sprintf('Correct selection %%'),'fontsize',fsy);
ylim([0 115]);
ytick = 0:20:100;
set(gca,'ytick',ytick); 

text(xsA,ysA,'D','fontsize',fsA,'Unit','normalized','fontname','Calibri');

lg = legend(methods_rfx,'location','northwest','fontsize',fsl,'Orientation','vertical');
set(lg,'linewidth',1);

%------------
subplot(nr,nc,5);
errorbarKxN(auc,auc*0,xgroups,'',colmap,0,bwe);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
title(sprintf('Model selection\nperformance'),'fontsize',fst,'fontname',fnt);
ylabel('Area under ROC','fontsize',fsy);
xlabel('Number of subjects','fontsize',fsy);
ylim([0.5 1.05]);ytick = .5:0.1:1;
set(gca,'ytick',ytick);
text(xsA,ysA,'E','fontsize',fsA,'Unit','normalized','fontname','Calibri');
%------------

subplot(nr,nc,6);
errorbarKxN(mdx,[eldx; ehdx],xgroups,'',colmap,0,bwe);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
title(sprintf('Parameter estimation\nperformance'),'fontsize',fst,'fontname',fnt);
ylabel('Error','fontsize',fsy);
xlabel('Number of subjects','fontsize',fsy);
ym = get(gca,'ylim'); ym = ym(2);
ylim([0 1.2*ym]); ytick = 0:.2:ym;
set(gca,'ytick',ytick);
text(xsA,ysA,'F','fontsize',fsA,'Unit','normalized','fontname','Calibri');

lg = legend(methods,'location','north','fontsize',fsl,'Orientation','horizontal');
set(lg,'linewidth',1);
end
