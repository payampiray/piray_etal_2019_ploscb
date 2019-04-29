function sim_B1(i)
% outliers in evidence space

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

run_out(nsim,1,ii);
run_out(nsim,2,ii);
run_out(nsim,3,ii);
end

function run_out(nsim,nout,ii)

simcat      = 'sim_B';

symname     = 'rwkf';
simstr      = sprintf('%s_outlier%d',symname,nout);

models      = {'model_rw1','model_kf1'};

pnames      = { {'\alpha','\beta'},{'\omega','\beta'} };
modelnames  = {'Reinforcement learning','Kalman filter'}; 
mnames      = {'RL','KF'}; 

normx1 = {@(x)1./(1+exp(-x)), @exp};
normx2 = {@exp,@exp};

normx  = {normx1,normx2};

sim_sim_B1(ii,simcat,simstr,models);
if ii(end)~=nsim, return; end
sim_wrap(nsim,simcat,simstr,modelnames,mnames,pnames,normx);
end

function sumfig
simcat     = 'sim_A';
fitsimstrs = {'fit_rwkf[10 30]'};
fsimfit0   = fullfile(getdefaults('pipedir'),simcat,fitsimstrs);

simcat     = 'sim_B';
fitsimstrs = {'fit_rwkf_outlier1','fit_rwkf_outlier2','fit_rwkf_outlier3'};
fsimfits   = [fsimfit0 fullfile(getdefaults('pipedir'),simcat,fitsimstrs)];

nout = 0:3;
[corrBMS,methods] = sim_stat_B1(fsimfits);

% close all;
fig_plot_evidence(corrBMS,nout,methods);

set(gcf,'name',mfilename);
end

function fig_plot_evidence(corrBMSs,xgroups,methods)
[fs,fst,fsl,fsy,fsalpha,fsxt,fsA,fn,fnt,xsA,ysA,fpos0,siz0,colmap,alf,cmaphbi,bw,colmapsim]=fig_plot_properties; %#ok<*ASGLU>

fpos = [.3 .35];
siz = fpos./fpos0(3:4).*siz0;

bw = .7*bw;

figure;
set(gcf,'units','centimeters');
fsiz = get(gcf,'position');
fsiz(3:4) = siz;
fsiz(1:2) = fpos0(1:2).*fsiz(1:2);
set(gcf,'position',fsiz);
set(gcf,'units','normalized');
%---------

mx = ceil(max(max(corrBMSs)));
ms = mx/5;

errorbarKxN(corrBMSs',0*corrBMSs',xgroups,'',colmap,0,bw);
set(gca,'fontsize',fs,'fontname',fn);
ylabel('Correct model selection %','fontsize',fsy);
title('Robustness of model selection to outliers','fontsize',fst,'fontname',fnt);
alpha(gca,alf);
set(gca,'ytick',0:ms:mx);
set(gca,'yticklabel',0:ms:mx);
ylim([0 1.2*mx]);
ax = ancestor(gca, 'axes');
xaxes = get(ax,'XAxis');
set(xaxes,'fontsize',fsxt);
legend(methods,'location','north','fontsize',fsl,'Orientation','horizontal');
xlabel('Number of outliers','fontsize',fsy); 

end