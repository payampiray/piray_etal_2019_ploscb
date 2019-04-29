function sim_C2(i)
% Figure 4
% Different ratio, HBI vs NHI

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
if nargin<1, ii=nsim; end

run_Nbar([20 20],nsim,ii);
run_Nbar([30 10],nsim,ii)
end

function run_Nbar(Nbar,nsim,ii)

simcat      = 'sim_C';

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
simcat   = 'sim_C';
fitsimstrs = {'fit_sim_rw[10 30]','fit_sim_rw[20 20]','fit_sim_rw[30 10]'}; 
xgroups = {'10/30','20/20','30/10'};
kref = 1;
ratek = [.25 .5 .75];

fsimfits  = fullfile(getdefaults('pipedir'),simcat,fitsimstrs);

[mpxp,spxp,mNbar,sNbar,auc,ms50,ms95,TR,ratings,xa,ya,auca,methods] = sim_stat_C2(fsimfits,kref);

modelrefname = 'RL';

mNbar =[mNbar; ratek];
sNbar =[sNbar; 0*ratek];
methods = [methods 'true'];

ratings = {'TA','FA','Accuracy'};

% close all;
fig_plot(modelrefname,xgroups,mpxp,spxp,mNbar,sNbar,ms50,ms95,TR,ratings,xa,ya,auca,methods);
set(gcf,'name',mfilename);
end

function fig_plot(modelrefname,xgroups,mpxp,spxp,mNbar,sNbar,ms50,ms95,TR,ratings,xa,ya,auca,methods)
[fs,fst,fsl,fsy,fsalpha,fsxt,fsA,fn,fnt,xsA,ysA,fpos0,siz0,colmap,alf,cmaphbi,bw,colmapsim]=fig_plot_properties; %#ok<*ASGLU>

fpos = [.4 .85];
siz = fpos./fpos0(3:4).*siz0;

xsA = -.2;
ysA = 1.15;

bwf = .9*bw;
bwi = 1.1*bw;
colmap = colmap([1 3],:);
colmapnbar = [colmap;colmapsim(2,:)];

nr = 3;
nc = 2;

figure;
set(gcf,'units','centimeters');
fsiz = get(gcf,'position');
fsiz(3:4) = siz;
fsiz(1:2) = fpos0(1:2).*fsiz(1:2);
set(gcf,'position',fsiz);
set(gcf,'units','normalized');
%---------

mpxp(mpxp<.01)=0.01;
ms95(ms95<1)=1;

subplot(nr,nc,1);
errorbarKxN(mpxp,spxp,xgroups,'',colmap,0,bw);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
title('Exceedance probability','fontsize',fst,'fontname',fnt);
ylim([0 1.2]);
ytick = 0:.2:1;
set(gca,'ytick',ytick);
ylabel(sprintf('PXP of %s',modelrefname),'fontsize',fsy);
% xlabel('#subjects per model','fontsize',fsy);
% legend(methods_rfx,'location','north','fontsize',fsl,'Orientation','horizontal');
% ax = ancestor(gca, 'axes');
% xaxes = get(ax,'XAxis');
% set(xaxes,'fontsize',fsxt);
text(xsA,ysA,'A','fontsize',fsA,'Unit','normalized','fontname','Calibri');

legend(methods(1:2),'location','north','fontsize',fsl,'Orientation','horizontal');
% ax = ancestor(gca, 'axes');
% xaxes = get(ax,'XAxis');
% set(xaxes,'fontsize',fsxt);
%-----------
subplot(nr,nc,2);
errorbarKxN(mNbar,sNbar,xgroups,'',colmapnbar,0,bwf);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
title(sprintf('Model comparison'),'fontsize',fst,'fontname',fnt);
% ylabel('Assigned subjects %','fontsize',fsy);
ylabel(sprintf('Frequency of %s',modelrefname),'fontsize',fsy);
% xlabel('#subjects per model','fontsize',fsy);
ylim([0 1.2]); ytick = 0:.2:1;
set(gca,'ytick',ytick);
% legend(methods,'location','north','fontsize',fsl,'Orientation','horizontal');
% ax = ancestor(gca, 'axes');
% xaxes = get(ax,'XAxis');
% set(xaxes,'fontsize',fsxt);
% text(xsA,ysA,'B','fontsize',fsA,'Unit','normalized','fontname','Calibri');

legend(methods,'location','north','fontsize',fsl,'Orientation','horizontal');
% ax = ancestor(gca, 'axes');
% xaxes = get(ax,'XAxis');
% set(xaxes,'fontsize',fsxt);
text(xsA,ysA,'B','fontsize',fsA,'Unit','normalized','fontname','Calibri');
%-----------
subplot(nr,nc,3);
errorbarKxN(ms50,ms50*0,xgroups,'',colmap,0,bw);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
title(sprintf('Model selection at PXP>0.5'),'fontsize',fst,'fontname',fnt);
ylabel(sprintf('Correct selection %%'),'fontsize',fsy);
xlabel('#subjects per model','fontsize',fsy);
ylim([0 110]);
text(xsA,ysA,'C','fontsize',fsA,'Unit','normalized','fontname','Calibri');
%-----------
subplot(nr,nc,4);
errorbarKxN(ms95,ms95*0,xgroups,'',colmap,0,bw);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
title(sprintf('Model selection at PXP>0.95'),'fontsize',fst,'fontname',fnt);
ylabel(sprintf('Correct selection %%'),'fontsize',fsy);
xlabel('#subjects per model','fontsize',fsy);
ylim([0 110]);
text(xsA,ysA,'D','fontsize',fsA,'Unit','normalized','fontname','Calibri');
%-----------
subplot(nr,nc,5);
errorbarKxN(TR,TR*0,ratings,'',colmap,0,bw);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
title(sprintf('Model attribution at r>0.95'),'fontsize',fst,'fontname',fnt);
ylabel(sprintf('Rate'),'fontsize',fsy);
% xlabel('#subjects per model','fontsize',fsy);
ylim([0 1.1]);
% ytick = 0:.2:1;
% set(gca,'ytick',ytick);
% legend(methods,'location','north','fontsize',fsl,'Orientation','horizontal');
ax = ancestor(gca, 'axes');
xaxes = get(ax,'XAxis');
set(xaxes,'fontsize',fsxt);
text(xsA,ysA,'E','fontsize',fsA,'Unit','normalized','fontname','Calibri');

%-----------
subplot(nr,nc,6);
for j=1:size(xa,1)
h1 = plot(xa{j},ya{j},'color',colmap(j,:),'linewidth',3); hold on;
h1.Color(4) = alf;
end    
set(gca,'fontsize',fs,'fontname',fn);
title(sprintf('Model attribution performance'),'fontsize',fst,'fontname',fnt);
xlabel('FA rate','fontsize',fsy);
ylabel('TA rate','fontsize',fsy);
set(gca,'ytick',0:0.2:1);
set(gca,'xtick',0:0.2:1);
% xlim([0 1.1]);
ylim([0 1.1]);
box off;

text(xsA,ysA,'F','fontsize',fsA,'Unit','normalized','fontname','Calibri');
% legend(methods,'location','southeast','fontsize',fsl,'Orientation','vertical');

basepos = get(gca,'position');
posloc = [1 1].*basepos(1:2) + [.7 .1].*(basepos(3:4));
posize = [.25 .5].*(basepos(3:4));
axes('Position',[posloc posize]);
errorbarKxN(auca,auca*0,{''},'',colmap,0,bwi);
ylim([.5 1]);
set(gca,'ytick',.5:.1:1);
alpha(gca,alf);
set(gca,'fontsize',fs,'fontname',fn);
% title(sprintf('MA performance'),'fontsize',fst,'fontname',fnt);
ylabel(sprintf('AUC'),'fontsize',fsy);
% lg = legend(methods,'location','north','fontsize',fsl,'Orientation','vertical');

end