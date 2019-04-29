function sim_C3(i)
% Supplementary Figure 1
% Non-hierarchical inference as a function of the variance of prior

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
vstr = {'V025','V1','V3','V6','V9'};
V    = [.25 1 3 6.25 9];

for i=1:5
    run_Nbar(ii,[10 30],V(i),vstr{i});
    run_Nbar(ii,[30 10],V(i),vstr{i});
end

end

function run_Nbar(ii,Nbar,v0,vstr)
nsim        = 100;
if nargin<1, ii=nsim; end

% % Nbar        = [10 30];

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

tx1   = [.1 1]';
tx2   = [.8 .4 3]';

logit = @(x)log(x./(1-x));

mu1  = [logit(tx1(1)); log(tx1(2))];
mu2  = [logit(tx2(1:2)); log(tx2(3))];
mu   = {mu1,mu2};

normx1 = {@(x)1./(1+exp(-x)), @exp};
normx2 = {@(x)1./(1+exp(-x)),@(x)1./(1+exp(-x)),@exp};
normx  = {normx1,normx2};

sim_sim(ii,simcat,simstr,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx);


m0          = {[0;0],[0;0;0]};
% v0          = 1; vstr = 'V1';
simstr0     = simstr;
simstr      = sprintf('%s%s',simstr0,vstr);
sim_sim_C3(ii,simcat,simstr,simcat,simstr0,m0,v0);

if ii(end)~=nsim, return; end
sim_wrap(nsim,simcat,simstr,modelnames,mnames,pnames,normx);

end

function sumfig
simcat   = 'sim_C';
fsum  = fullfile(getdefaults('sumdir'),simcat,sprintf('sum_C3.mat'));
if ~exist(fsum,'file')
    simstrs(1,:)  = {'fit_sim_rw[10 30]V025','fit_sim_rw[10 30]V1','fit_sim_rw[10 30]V3','fit_sim_rw[10 30]V625','fit_sim_rw[10 30]V9'};
    simstrs(2,:)  = {'fit_sim_rw[30 10]V025','fit_sim_rw[30 10]V1','fit_sim_rw[30 10]V3','fit_sim_rw[30 10]V625','fit_sim_rw[30 10]V9'};
    xgroups = {'0.25', '1', '3', '6.25', '9'};    
    fsimfits  = fullfile(getdefaults('pipedir'),simcat,simstrs );        

    [mx,ex,mp,ep,methods] = sim_stat_C3(fsimfits);    
    save(fsum,'xgroups','mx','ex','mp','ep')
else
    sums    = load(fsum);
    mx      = sums.mx;
    ex      = sums.ex;
    mp      = sums.mp;
    ep      = sums.ep;
    xgroups = sums.xgroups;
end

% close all;
fig_plot(xgroups,mx,ex,mp,ep);
set(gcf,'name',mfilename);
end

function fig_plot(xgroups,mx,ex,mp,ep)

[fs,fst,fsl,fsy,fsalpha,fsxt,fsA,fn,fnt,xsA,ysA,fpos0,siz0,colmap,alf,cmaphbi,bw,colmapsim]=fig_plot_properties; %#ok<*ASGLU>
cmaplap = colmap(1,:);

fpos = [.5 fpos0(4)];
siz = fpos./fpos0(3:4).*siz0;

bwp = 1.4*bw;

ysA = ysA+.05;

nr = 2;
nc = 2;

figure;
set(gcf,'units','centimeters');
fsiz = get(gcf,'position');
fsiz(3:4) = siz;
fsiz(1:2) = fpos0(1:2).*fsiz(1:2);
set(gcf,'position',fsiz);
set(gcf,'units','normalized');
%---------
abc = 'ABCD';

for i=1:2
subplot(nr,nc,(i-1)+1 );
errorbarKxN(mp{i},ep{i},'','',cmaplap,0,bwp);
set(gca,'fontsize',fs,'fontname',fn);
set(gca,'xticklabel',xgroups,'xcolor','k');
% set(gca,'ygrid','on');
alpha(gca,alf);
ylabel(sprintf('PXP'),'fontsize',fsy);
ylim([0 1.05]);
set(gca,'ytick',0:0.2:1);
title(sprintf('Scenario %d\n Model selection of NHI',i),'fontsize',fst,'fontname',fnt);
ax = ancestor(gca, 'axes');
xaxes = get(ax,'XAxis');
set(xaxes,'fontsize',fsxt);
text(xsA,ysA,abc(i-1+1 ) ,'fontsize',fsA,'Unit','normalized','fontname','Calibri');

subplot(nr,nc,(i-1)+3 );
errorbarKxN(mx{i},ex{i},'','',cmaplap,0,bwp);
set(gca,'fontsize',fs,'fontname',fn);
set(gca,'xticklabel',xgroups,'xcolor','k');
% set(gca,'ygrid','on');
alpha(gca,alf);
ylabel(sprintf('Error'),'fontsize',fsy);
% ylim([0 1.05]);
% set(gca,'ytick',0:0.2:1);
ylabel(sprintf('Error of %s','\alpha'),'fontsize',fst,'fontname',fnt);
title(sprintf('Parameter estimation of NHI'),'fontsize',fst,'fontname',fnt);
ax = ancestor(gca, 'axes');
xaxes = get(ax,'XAxis');
set(xaxes,'fontsize',fsxt);

xlabel('Prior variance');

text(xsA,ysA,abc(i-1+3 ) ,'fontsize',fsA,'Unit','normalized','fontname','Calibri');
end
%-----

end
