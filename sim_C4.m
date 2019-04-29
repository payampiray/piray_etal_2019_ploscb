function sim_C4(i)
% Supplementary Figure 2
% The same as sim_C1, but with a different learning rate setting

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
nsim        = 20;

simcat      = 'sim_C';
Nbar        = [10 30];

symname     = 'sim_rw2';
simstr      = sprintf('%s%s',symname,mat2str(Nbar));


models = {'model_rw1','model_rw2'};
tasks  = {'task_rw1','task_rw2'};
pnames = { {'\alpha','\beta'},{'\alpha^+','\alpha^-','\beta'} };

modelnames = {'Regular RL','Dual-\alpha RL'}; 
mnames = modelnames;

tau1 = 4*[1 1]';
tau2 = 4*[1 1 1]';
tau  = {tau1,tau2};

% % rw2
% tx1  = [.6 1]';
% tx2  = [.8 .4 3]'; 
% rw3
tx1  = [.6 1]';
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
simstr   = 'sim_rw2[10 30]';
fsimfit  = fullfile(getdefaults('pipedir'),simcat,sprintf('fit_%s.mat',simstr));

simfit   = load(fsimfit);
fit      = simfit.fit;
config   = simfit.config;
mnames     = config.mnames;
pnames     = config.pnames;

mnames{1} = 'RL';

[pxpm,pxps,Nbarm,Nbars,ngm,egm,rmm,rms,mdx,sdx,corrBMS,methods] = sim_stat_A1(fit,config);

Nbarm =[Nbarm; .25 .75];
Nbars =[Nbars; 0 0];
Nbarstrs = {'HBI','True'};

% close all;
fig_plot(pxpm,pxps,Nbarm,Nbars,ngm,egm,rmm,rms,mdx,sdx,corrBMS,Nbarstrs,mnames,pnames,methods);
set(gcf,'name',mfilename);
end

function fig_plot(pxpm,pxps,Nbar,Nbars,mgm,egm,rm,re,mdx,sdx,corrBMS,Nbarstrs,mnames,pnames,methods)
pnames = sim_plot_adjust_pnames(pnames);

[fs,fst,fsl,fsy,fsalpha,fsxt,fsA,fn,fnt,xsA,ysA,fpos0,siz0,colmap,alf,cmaphbi,bw,colmapsim]=fig_plot_properties; %#ok<*ASGLU>
cmaphbi2 = [cmaphbi; colmapsim(2,:)];

fpos = [.6 fpos0(4)];
siz = fpos./fpos0(3:4).*siz0;

xst = .5;
yst = 1.1;
ysr = .1;

bwp = .9*bw;
bwr = .9*bw;
bwi = 1.6*bw;
bws = 1.4*bw;
bwe = [bw*.5 bw*.8];

nr = 2;
nc = 3;

figure;
set(gcf,'units','centimeters');
fsiz = get(gcf,'position');
fsiz(3:4) = siz;
fsiz(1:2) = fpos0(1:2).*fsiz(1:2);
set(gcf,'position',fsiz);
set(gcf,'units','normalized');
%---------


pxpm(pxpm<0.01)=0.01;
h(1) = subplot(nr,nc,1);
errorbarKxN(pxpm,pxps,mnames,'',cmaphbi,0,bwp);
set(gca,'fontsize',fs,'fontname',fn);
set(gca,'xticklabel',mnames,'xcolor','k');
% set(gca,'ygrid','on');
alpha(gca,alf);
ylabel(sprintf('PXP'),'fontsize',fsy);
ylim([0 1.05]);
set(gca,'ytick',0:0.2:1);
title(sprintf('Model selection by HBI'),'fontsize',fst,'fontname',fnt);
ax = ancestor(gca, 'axes');
xaxes = get(ax,'XAxis');
set(xaxes,'fontsize',fsxt);
text(xsA,ysA,'A','fontsize',fsA,'Unit','normalized','fontname','Calibri');
    
% figure;
h(2) = subplot(nr,nc,2);
errorbarKxN(Nbar,Nbars,mnames,'',cmaphbi2,0,bwp);
set(gca,'fontsize',fs,'fontname',fn);
set(gca,'xticklabel',mnames,'xcolor','k');
% set(gca,'ygrid','on');
alpha(gca,alf);
ylim([0 1.05])
ylabel(('Model frequency'),'fontsize',fsy);
% title(sprintf('Estimated frequency of models\n across subjects'),'fontsize',fst,'fontname',fnt);
title(sprintf('Model comparison by HBI'),'fontsize',fst,'fontname',fnt);
ax = ancestor(gca, 'axes');
xaxes = get(ax,'XAxis');
set(xaxes,'fontsize',fsxt);
text(xsA,ysA,'B','fontsize',fsA,'Unit','normalized','fontname','Calibri');
lg = legend(Nbarstrs,'location','northwest','fontsize',fsl,'Orientation','vertical'); hold on;
set(lg,'linewidth',1);


h(3) = subplot(nr,nc,3);
errorbarKxN(rm,re,{'TA','FA'},'',cmaphbi,0,bwr);
set(gca,'fontsize',fs,'fontname',fn);
% set(gca,'xticklabel',{'Correctly assigned','Incorrectly assigned'},'xcolor','k');
% set(gca,'xticklabel',{'True Positive','False Positive'},'xcolor','k','XTickLabelRotation',20);
ylabel('Responsibility','fontsize',fsy);
title('Model attribution by HBI','fontsize',fst,'fontname',fnt);
ax = ancestor(gca, 'axes');
xaxes = get(ax,'XAxis');
set(xaxes,'fontsize',fsxt);
alpha(gca,alf);
% set(gca,'ygrid','on');
set(gca,'ytick',.5:0.1:1);
% set(gca,'yticklabel',.5:.1:1);
ylim([.5 1.025]);
text(xsA,ysA,'C','fontsize',fsA,'Unit','normalized','fontname','Calibri');


%------------------
basepos = get(gca,'position');
pos = [1 1].*basepos(1:2) + [.5 .75].*(basepos(3:4));
hinset = axes('Position',[pos .1 0.06]);
bar(mgm*100,'Horizontal','on','facecolor',cmaphbi,'EdgeColor','k','linewidth',1,'basevalue',50,'barwidth',bwi); hold on;
plot([mgm-egm mgm+egm]*100,[1 1],'-','color','k','linewidth',2);
xlim([50 100]);
% set(gca,'box','off');
set(hinset,'xtick',50:10:90);
set(hinset,'xticklabel',50:10:90);
set(hinset,'ticklength', [0 0]);
set(hinset,'yticklabel',[],'xcolor','k');
set(hinset,'linewidth',1);
title('Accuracy %','fontsize',fsxt,'fontname',fnt);
alpha(gca,alf);
%--------------------

%--
% colmap = [.2 .2 .2; .5 .5 .5; .9 .9 .9];

% plotbarit(corrBMS,0*corrBMS,methods);
h(4) = subplot(2,3,4);
corrBMS = max(corrBMS,1);
errorbarNxK(corrBMS,corrBMS*0,methods,'',colmap,0,bws);
set(gca,'fontsize',fs,'fontname',fn);
ylabel('Correct selection %','fontsize',fsy);
% title(sprintf('Model selection\nperformance'),'fontsize',fst,'fontname',fnt);
set(gca,'ytick',0:20:100);
ylim([0 110]);
ax = ancestor(gca, 'axes');
xaxes = get(ax,'XAxis');
set(xaxes,'fontsize',fsxt);
text(xsA,ysA+ysr,'D','fontsize',fsA,'Unit','normalized','fontname','Calibri');
alpha(gca,alf);

text(xst,yst,sprintf('Model selection\nperformance'),'fontsize',fst,'fontname',fnt,...
    'Unit','normalized','fontweight','bold','Parent',h(4),'HorizontalAlignment','center'); hold on;   

% figure for parameters

% figure for parameters
EF ={'E','F'};
dls = [.05 .04];
if ~isempty(mdx)    
for k=[1 2]
    hk = subplot(2,3,4+k);

%     h(4+ki) = errorbarKxN(mdx{k},sdx{k},pnames{k},'',colmap);
%     errorbarKxN_multiaxes(mdx{k},sdx{k},[1 2],pnames{k},'',colmap);
%     set(gca,'fontsize',fs,'fontname',fn);
%     set(gca,'xticklabel',pnames{k},'xcolor','k','FontSmoothing','on');
    
    %----------------------------------------------------------------------
    pos = get(hk,'Position');
    text(xst,yst,sprintf('Parameter estimation\n performance for %s',mnames{k}),'fontsize',fst,'fontname',fnt,...
        'Unit','normalized','fontweight','bold','Parent',hk,'HorizontalAlignment','center'); hold on;    
    set(hk,'visible','off');
    text(xsA,ysA+ysr,EF{k},'fontsize',fsA,'Unit','normalized','fontname','Calibri','Parent',hk); hold on;
    
    np = size(mdx{k},2);
    dl = dls(k);
    
    pos1 = pos;
    x0 = pos1(1);    
    lp = 1/np* (pos1(3) -(np-1)*dl);

    for i=1:np
        pos1 = pos;
        pos1(1) = x0+(i-1)*(lp+dl);
        pos1(3) = lp;

        axes('Position',pos1);    
%         errorbarKxN(mxf(:,i),exf(:,i),pnamesf(i),'',cmaphbi,0,bw*rbf); hold on;
        errorbarKxN(mdx{k}(:,i),sdx{k}(:,i),pnames{k}{i},'',colmap,0,bwe(k)); hold on;        
        alpha(gca,alf);
        set(gca,'fontsize',fs,'fontname',fn);

        ax = ancestor(gca, 'axes');
        xaxes = get(ax,'XAxis');
        set(xaxes,'Fontsize',fsalpha,'TickLabelInterpreter','latex');
        if ~strcmp(pnames{k}{i}(1:3),'\it')
            set(xaxes,'FontWeight','bold');
        end
        
        if i==1
            ylabel('Error','fontsize',fsy);
        end

        ytick = get(gca,'ytick');
        ytick(end) = [];
        ytick = ytick(1:2:end);
        set(gca,'ytick',ytick);


    end    
end

end

yl = get(gca,'ylim');
yl(2) = 1.4*yl(2);
ylim(yl);
lg = legend(methods,'location','north','fontsize',fsl,'Orientation','vertical'); hold on;
set(lg,'linewidth',1);
poslg = get(lg,'position');
% poslg(1) = poslg(1)+.04;
set(lg,'position',poslg);
end
