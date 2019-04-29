function sim_H3(i)
% Figure 12
% HBI t-test for skewed distribution

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

skewness    = -.5;
simcat      = 'sim_H3';
Nbar        = 50;

symname     = 'brl';
simstr      = sprintf('%s%d[%+0.2f]',symname,Nbar,skewness);


models     = {'model_rwneut1'};
tasks      = {'task_rwgo1'};
modelnames = {'Biased RL'};
mnames     = {'Biased RL'};
pnames     = { {'\alpha','\beta','\it b'}};

tau1 = 4*[1 1 nan]';
tau  = {tau1};

tx1  = [.1 2 nan]';
logit = @(x)log(x./(1-x));

mu1  = [logit(tx1(1)); log(tx1(2)); tx1(3)];
mu   = {mu1};

fg1  = [1;1;0];
fg   = {fg1};

normx1 = {@(x)1./(1+exp(-x)),@exp};
normx  = {normx1};

% note: kur=3 is the kur of normal (the 4th input in persrnd)
fanorm = @(n)pearsrnd(0,1,skewness,3,1,n);

anormal1 = {[],[],fanorm };

anormal = {anormal1};

sim_sim_H(ii,simcat,simstr,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx,fg,anormal);
if ii(end)~=nsim, return; end
sim_wrap(nsim,simcat,simstr,modelnames,mnames,pnames,normx);
end

function sumfig
simcat = 'sim_H3';
fitsimstrs(1,:) = {'fit_brl20[-0.50].mat','fit_brl50[-0.50].mat'}; 
NN = [20 50];

mref = 0;
fsimfits  = fullfile(getdefaults('pipedir'),simcat,fitsimstrs);

[tr,xb,xl,xh,xc,sighbi,sigbarhbi,sigsim,sigbarsim] = sim_stat_H3(fsimfits,mref);

nh = 10^7;
[ssamples] = pearsrnd(0,1,-.5,3,1,nh);

% close all;
plot_fig(ssamples,tr,xb,xl,xh,sighbi,sigbarhbi,sigsim,sigbarsim,NN);
set(gcf,'name',mfilename);
end

function plot_fig(ssamples,tr,xb,xl,xh,sighbi,sigbarhbi,sigsim,sigbarsim,NN)

[fs,fst,fsl,fsy,fsalpha,fsxt,fsA,fn,fnt,xsA,ysA,fpos0,siz0,colmap,alf,cmaphbi,bw,colmapsim]=fig_plot_properties; %#ok<*ASGLU>

bw = .35*bw;

fpos = [.33 .84];
siz = fpos./fpos0(3:4).*siz0;

colmapsim = [cmaphbi; colmapsim(2,:)];

%--------

methods = {'HBI','Benchmark'};

figure;
set(gcf,'units','centimeters');
fsiz = get(gcf,'position');
fsiz(3:4) = siz;
fsiz(1:2) = fpos0(1:2).*fsiz(1:2);
set(gcf,'position',fsiz);
set(gcf,'units','normalized');

nr = 3;
nc = 2;

abc = 'ABCDEF';

ymax = 1.05*max(max(sigbarhbi));

subplot(nr,nc,1); 
pos = get(gca,'position');

hk = subplot(nr,nc,1:2); 
get(hk,'position');

hh = histogram(ssamples);   
set(hh,'FaceColor','k');
xlim([-5 5]);
set(gca,'ytick',[]);
set(gca,'fontsize',fs,'fontname',fn);
title('Skewed distribution','fontsize',fst,'fontname',fnt);

hk = axes('Position',pos);
set(hk,'visible','off');
text(xsA,ysA,abc(1),'fontsize',fsA,'Unit','normalized','fontname','Calibri','parent',hk);

%--
for j=1:2
    mx = [sighbi(j);sigsim(j)];
    
    exl = [sigbarhbi(1,j);sigbarsim(1,j)];
    exh = [sigbarhbi(2,j);sigbarsim(2,j)];
    ex  = [exl; exh];

    subplot(nr,nc,2+j);
    errorbarKxN(mx,ex,{sprintf('N=%d',NN(j))},'',colmapsim,0,bw); hold on;
    alpha(gca,alf);
    set(gca,'fontsize',fs,'fontname',fn);

    ax = ancestor(gca, 'axes');
    xaxes = get(ax,'XAxis');
    set(xaxes,'Fontsize',fsalpha);

    ylim([0 ymax]);
    ytick = get(gca,'ytick');
    set(gca,'ytick',ytick);
    
    title('Inference at P<0.05','fontsize',fst,'fontname',fnt);

    if j==1
        hyl =  ylabel(sprintf('Probability'),'fontsize',fsy);
    end        
    text(xsA,ysA,abc(1+j),'fontsize',fsA,'Unit','normalized','fontname','Calibri');        
end

    lg = legend(methods,'location','north','fontsize',fsl,'Orientation','vertical'); 
    set(lg,'linewidth',1);    

%------------------
for i=1:2
    subplot(nr,nc, 4+i);    
    xbj = xb;
    xlj = xl;
    xhj = xh;

    tt = 2:size(xbj,2)-1;
    plot(tr,xbj(i,:),'color',cmaphbi,'linewidth',2); hold on;
    plot(0:.01:1,0:.01:1,'k','linewidth',2); hold on;
    errorbar(tr(tt),xbj(i,tt),xbj(i,tt)-xlj(i,tt),xhj(i,tt)-xbj(i,tt),'color',cmaphbi,'linewidth',2); hold on;
    
    set(gca,'fontsize',fs,'fontname',fn);
    title(sprintf('Performance for N=%d',NN(i)),'fontsize',fst,'fontname',fnt);
    
    ylim([0 1]);
    ytick = 0:.2:.8;
    set(gca,'ytick',ytick);
    xlim([0 1]);
    xtick = 0:.2:1;
    set(gca,'xtick',xtick);
    
    xlabel('P-value','fontsize',fsy);
    
    if i==1    
    hylcdf = ylabel(sprintf('Probability\nunder the null'),'fontsize',fsy);
    set(hylcdf,'Units','centimeters')
    poshyl = get(hylcdf,'position');
    
    set(hyl(1),'Units','centimeters')
    poshyl1 = get(hyl(1),'position');
    poshyl(1) = poshyl1(1);
    set(hylcdf,'position',poshyl);
        
    end
    if i==2
    lg = legend({'HBI','Theory'},'location','southeast','fontsize',fsl,'Orientation','vertical'); 
    set(lg,'linewidth',1);    
    end
    text(xsA,ysA,abc(3+i),'fontsize',fsA,'Unit','normalized','fontname','Calibri');    
end

end
