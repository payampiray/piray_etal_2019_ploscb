function sim_F1(i)
% Figure 7
% outliers in parameters space

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

run0_N(30,ii,nsim);
for n=1:4
    Nbar = [30 n];
    run1_N(Nbar,ii,nsim);
    run2_N(Nbar,ii,nsim);
    run3_N(Nbar,ii,nsim);
end

end

function run0_N(Nbar,ii,nsim)

simcat      = 'sim_F'; % previously outsimple2

symname     = 'rl0';
simstr      = sprintf('%s%s',symname,mat2str([Nbar 0]));


models = {'model_rw1'};
tasks  = {'task_rw1'};
pnames = { {'\alpha','\beta'} };

modelnames = {'Regular RL'};
mnames = modelnames;

tau1 = 4*[1 1]';
tau  = {tau1};

tx1  = [.3 3]';
logit = @(x)log(x./(1-x));

mu1  = [logit(tx1(1)); log(tx1(2))];

mu   = {mu1};

normx1 = {@(x)1./(1+exp(-x)), @exp};

normx  = {normx1};


sim_sim_F1(ii,simcat,simstr,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx);
if ii(end)~=nsim, return; end
sim_wrap(nsim,simcat,simstr,modelnames,mnames,pnames,normx);
end

function run1_N(Nbar,ii,nsim)

simcat      = 'sim_F'; % previously outsimple2

symname     = 'rl1';
simstr      = sprintf('%s%s',symname,mat2str(Nbar));


models = {'model_rw1','model_rw1'};
tasks  = {'task_rw1','task_rw1'};
pnames = { {'\alpha','\beta'},{'\beta'} };

modelnames = {'Regular RL','Null'};
mnames = modelnames;

tau1 = 4*[1 1]';
tau2 = [inf inf]';
tau  = {tau1,tau2};

tx1  = [.3 3]';
tx2  = [.99 1];
logit = @(x)log(x./(1-x));

mu1  = [logit(tx1(1)); log(tx1(2))];
mu2  = [logit(tx2(1)); log(tx2(2))];

mu   = {mu1,mu2};

normx1 = {@(x)1./(1+exp(-x)), @exp};
normx2 = {@(x)1./(1+exp(-x)), @exp};

normx  = {normx1,normx2};


sim_sim_F1(ii,simcat,simstr,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx);
if ii(end)~=nsim, return; end
sim_wrap(nsim,simcat,simstr,modelnames,mnames,pnames,normx);
end

function run2_N(Nbar,ii,nsim)

simcat      = 'sim_F'; % previously outsimple2

symname     = 'rl2';
simstr      = sprintf('%s%s',symname,mat2str(Nbar));


models = {'model_rw1','model_rw1'};
tasks  = {'task_rw1','task_rw1'};
%%-
pnames = { {'\alpha','\beta'},{'\alpha','\beta'} };

modelnames = {'Regular RL','Regular RL'}; 
mnames = modelnames;

tau1 = 4*[1 1]';
tau2 = inf*[1 1]';
tau  = {tau1,tau2};

tx1  = [.3 3]';
tx2  = [.3 .1]';
logit = @(x)log(x./(1-x));

mu1  = [logit(tx1(1)); log(tx1(2))];
mu2  = [logit(tx2(1)); log(tx2(2))];

mu   = {mu1,mu2};
%%-

normx1 = {@(x)1./(1+exp(-x)), @exp};
normx2 = {@exp};

normx  = {normx1,normx2};


sim_sim_F1(ii,simcat,simstr,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx);
if ii(end)~=nsim, return; end
sim_wrap(nsim,simcat,simstr,modelnames,mnames,pnames,normx);
end

function run3_N(Nbar,ii,nsim)

simcat      = 'sim_F';

symname     = 'rl3';
simstr      = sprintf('%s%s',symname,mat2str(Nbar));


models = {'model_rw1','model_rw1'};
tasks  = {'task_rw1','task_rw1'};
%%-
pnames = { {'\alpha','\beta'},{'\alpha','\beta'} };

modelnames = {'Regular RL','Regular RL'};
mnames = modelnames;

tau1 = 4*[1 1]';
tau2 = inf*[1 1]';
tau  = {tau1,tau2};

tx1  = [.3 3]';
tx2  = [.1 .1]';
logit = @(x)log(x./(1-x));

mu1  = [logit(tx1(1)); log(tx1(2))];
mu2  = [logit(tx2(1)); log(tx2(2))];

mu   = {mu1,mu2};
%%-

normx1 = {@(x)1./(1+exp(-x)), @exp};
normx2 = {@exp};

normx  = {normx1,normx2};


sim_sim_F1(ii,simcat,simstr,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx);
if ii(end)~=nsim, return; end
sim_wrap(nsim,simcat,simstr,modelnames,mnames,pnames,normx);
end

function sumfig
simcat   = 'sim_F';
fsum    = fullfile(getdefaults('sumdir'),simcat,'sum.mat');
if ~exist(fsum,'file')
    fitsimstrs(1,:) = {'fit_rl0[30 0]','fit_rl2[30 1]','fit_rl2[30 2]','fit_rl2[30 3]','fit_rl2[30 4]'};
    fitsimstrs(2,:) = {'fit_rl0[30 0]','fit_rl3[30 1]','fit_rl3[30 2]','fit_rl3[30 3]','fit_rl3[30 4]'};
    xgroups = 0:4;
    kref = 1;

    fsimfits  = fullfile(getdefaults('pipedir'),simcat,fitsimstrs);
    [mdgx,sdgx,pnames,methods]=sim_stat_F1(fsimfits,kref);
    save(fsum,'xgroups','mdgx','sdgx','pnames','methods');
else
    sums    = load(fsum);
    xgroups = sums.xgroups;
    mdgx    = sums.mdgx;
    sdgx    = sums.sdgx;
    pnames  = sums.pnames;
    methods = sums.methods;
end

% close all;
fig(mdgx,sdgx,pnames, methods,xgroups);
set(gcf,'name',mfilename);
end

function fig(mdx,sdx,pnames, methods,xgroups)
[fs,fst,fsl,fsy,fsalpha,fsxt,fsA,fn,fnt,xsA,ysA,fpos0,siz0,colmap,alf,cmaphbi,bw,colmapsim]=fig_plot_properties; %#ok<*ASGLU>

fpos = [0.4 .3];
siz = fpos./fpos0(3:4).*siz0;

xsA = -.2;
ysA = 1.15;

figure;
set(gcf,'units','centimeters');
fsiz = get(gcf,'position');
fsiz(3:4) = siz;
fsiz(1:2) = fpos0(1:2).*fsiz(1:2);
set(gcf,'position',fsiz);
set(gcf,'units','normalized');
%---------

[nc,nr] = size(mdx);

ABC = {'A','B','C','D','E'};

yl = nan(nr,nc);

ylpos = nan(nr,3);

ki = 1;
for i=1:nr
    for j=1:nc
    hp(i,j) = subplot(nr,nc,ki);      %#ok<AGROW>
    
    for m=1:3
        hep = plot(xgroups,mdx{j,i}(m,:),'color',colmap(m,:),'linewidth',2); hold on;
        hep.Color(4) = alf;      
        for k=1:length(xgroups)
            ax = xgroups(k);
            el = mdx{j,i}(m,k) - sdx{j,i}(m,k);
            eh = mdx{j,i}(m,k) + sdx{j,i}(m,k);
            hep = plot([ax;ax],[el;eh],'-','color',colmap(m,:),'linewidth',2);
            hep.Color(4) = alf;
            
            axeps = ax + [-0.05 .05];
            plot(axeps,[el;el],'-','color',colmap(m,:),'linewidth',2);            
            plot(axeps,[eh;eh],'-','color',colmap(m,:),'linewidth',2);            
        end
    end    
    set(gca,'fontsize',fs,'fontname',fn);
    set(gca,'xtick',xgroups);
        
    xlim([-.1+xgroups(1) xgroups(end)+.1]);
    
    if i==1        
        text(xsA,ysA,ABC{j},'fontsize',fsA,'Unit','normalized','fontname','Calibri');
        title(sprintf('Scenario %d',j),'fontsize',fst,'fontname',fnt);
    end
    if i==2
        xlabel('Number of outliers','fontsize',fsy);
    end
    if j==1
        ylabel(sprintf('Error of %s',pnames{i}),'fontsize',fsy,'fontname',fnt,'interpreter','tex');        
        ylh = ylabel( ['Error of\fontsize{18} ',(pnames{i})],'fontsize',fsy,'fontname',fnt,'interpreter','tex');
        ylpos(i,:) = get(ylh,'position');
        if i==2
            ylpos(i,1)=ylpos(1,1);
            set(ylh,'position',ylpos(i,:));
        end
    end
    
    yll = get(gca,'ylim');
    yl(i,j) = 1.1*yll(2);

    ki=ki+1;
    end
end

for i=1:nr
    yll = [0 .99*max(yl(i,:))];
    for j=1:nc
        set(hp(i,j),'ylim',yll);
    end
end

hx = nan(1,3);
i=1; j=nc;
subplot(nr,nc,j);
for m=1:3
hx(m) = plot(xgroups,mdx{j,i}(m,:),'color',colmap(m,:),'linewidth',2);
colm = [colmap(m,:), alf];
set(hx(m),'Color',colm);
end
legend(hx,methods,'Orientation','horizontal','location','Northeast','fontsize',fsl);

end