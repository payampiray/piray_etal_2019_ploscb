function sim_H2(i)
% Figure 11
% HBI t-test under the null

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
nsim        = 2000;
if nargin<1, ii=nsim; end

run_b1(ii,nsim);
run_b2(ii,nsim);
end

function run_b1(ii,nsim)

simcat      = 'sim_H2';

symname     = 'brl';
simstr      = sprintf('%s[null]',symname);

Nbar   = 20;

models     = {'model_rwneut1'};
tasks      = {'task_rwgo1'};
modelnames = {'Biased RL'};
mnames     = {'Biased RL'};
pnames     = { {'\alpha','\beta','\it b'}};

tau1 = 4*[1 1 .25]';
tau  = {tau1};

tx1  = [.1 2 0]';
logit = @(x)log(x./(1-x));

mu1  = [logit(tx1(1)); log(tx1(2)); tx1(3)];
mu   = {mu1};

fg1  = [1;1;0];
fg   = {fg1};

normx1 = {@(x)1./(1+exp(-x)),@exp,@(x)x};
normx  = {normx1};

sim_sim_H(ii,simcat,simstr,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx,fg);
if ii(end)~=nsim, return; end
sim_wrap(nsim,simcat,simstr,modelnames,mnames,pnames,normx);
end

function run_b2(ii,nsim)

simcat      = 'sim_H2';

symname     = 'brl2';
simstr      = sprintf('%s[null]',symname);

Nbar   = [20 20];

models      = {'model_rwneut1','model_rw2'};
tasks       = {'task_rwgo1','task_rwgo2'};
modelnames  = {'Biased RL','Dual-\alpha RL'};
mnames      = modelnames;
pnames      = { {'\alpha','\beta','\it b'},{'\alpha^+','\alpha^-','\beta'} };

tau1 = 4*[1 1 .25]';
tau2 = 4*[1 1 1]';
tau  = {tau1,tau2};

tx1  = [.1 1 0]';
tx2  = [.8 .4 3]';
logit = @(x)log(x./(1-x));

mu1  = [logit(tx1(1)); log(tx1(2)); (tx1(3))];
mu2  = [logit(tx2(1:2)); log(tx2(3));];

mu   = {mu1,mu2};

fg1  = [1;1;0];
fg2  = [1;1;1];

fg   = {fg1,fg2};

normx1 = {@(x)1./(1+exp(-x)),@exp,@(x)x};
normx2 = {@(x)1./(1+exp(-x)),@(x)1./(1+exp(-x)),@exp};
normx  = {normx1;normx2};

sim_sim_H(ii,simcat,simstr,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx,fg);
if ii(end)~=nsim, return; end
sim_wrap(nsim,simcat,simstr,modelnames,mnames,pnames,normx);
end

function sumfig
simcat = 'sim_H2';
fitsimstrs(1,:) = {'fit_brl[null].mat','fit_brl2[null].mat'}; 
          
fsimfits  = fullfile(getdefaults('pipedir'),simcat,fitsimstrs);
[tr,xb,xl,xh,xc] = sim_stat_H2(fsimfits,0);

% close all;
fig(tr,xb,xl,xh);
set(gcf,'name',mfilename);
end

function fig(tr,xb,xl,xh)

[fs,fst,fsl,fsy,fsalpha,fsxt,fsA,fn,fnt,xsA,ysA,fpos0,siz0,colmap,alf,cmaphbi,bw,colmapsim]=fig_plot_properties; %#ok<*ASGLU>

% ysA = 1.13;
fpos = [.4 .33];
siz = fpos./fpos0(3:4).*siz0;

%--------
posm = 250; % pixel

figure;
set(gcf,'units','centimeters');
fsiz = get(gcf,'position');
fsiz(3:4) = siz;
fsiz(1:2) = fpos0(1:2).*fsiz(1:2);
set(gcf,'position',fsiz);
set(gcf,'units','normalized');

nr = 1;
nc = size(xb,1);

abc = 'AB';
ki = 1;
for i=1:nc
    subplot(nr,nc,ki);
    set(gca,'units','pixels');
    pos  = get(gca,'position');
    pos(3:4) = posm;
    set(gca,'position',pos);

    set(gca,'units','normalized');
    pos  = get(gca,'position');    
    pos(2) = pos(2)*1.3;
    set(gca,'position',pos);
    
    xbj = xb;
    xlj = xl;
    xhj = xh;

    tt = 2:size(xbj,2)-1;    
    plot(tr,xbj(i,:),'color',cmaphbi,'linewidth',2); hold on;
    plot(0:.01:1,0:.01:1,'k','linewidth',1);    
    errorbar(tr(tt),xbj(i,tt),xbj(i,tt)-xlj(i,tt),xhj(i,tt)-xbj(i,tt),'color',cmaphbi,'linewidth',2); hold on;

    set(gca,'fontsize',fs,'fontname',fn);
    title(sprintf('Performance in scenario %d',i),'fontsize',fst,'fontname',fnt);
    
    ylim([0 1]);
    ytick = 0:.2:1;
    set(gca,'ytick',ytick);
    xlim([0 1]);
    xtick = 0:.2:1;
    set(gca,'xtick',xtick);
    
    xlabel('P-value','fontsize',fsy);
    ylabel('Probability under the null','fontsize',fsy);
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname','Calibri');
    
    lg = legend({'HBI','Theory'},'location','southeast','fontsize',fsl,'Orientation','vertical');
    set(lg,'linewidth',1);
    
    ki = ki+1;
end
end