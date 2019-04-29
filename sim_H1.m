function sim_H1(i)
% Figure 10
% HBI t-test performance given different effect sizes

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

run_b1(+.25,ii,nsim);
run_b1(+.5 ,ii,nsim);
run_b1(+.75,ii,nsim);
run_b1(+1  ,ii,nsim);
% 
run_b2(+.25,ii,nsim);
run_b2(+.5 ,ii,nsim);
run_b2(+.75,ii,nsim);
run_b2(+1  ,ii,nsim);
end

function run_b1(b,ii,nsim)

simcat      = 'sim_H1';

symname     = 'brl';
simstr      = sprintf('%s[%+0.2f]',symname,b);

Nbar        = 20;

models      = {'model_rwneut1'};
tasks       = {'task_rwgo1'};
modelnames  = {'Biased RL'};
mnames      = {'Biased RL'};
pnames      = { {'\alpha','\beta','\it b'}};

tau1 = 4*[1 1 .25]';
tau  = {tau1};

tx1   = [.1 2 b]';
logit = @(x)log(x./(1-x));

mu1  = [logit(tx1(1)); log(tx1(2)); tx1(3)];
mu   = {mu1};

fg1  = [1;1;0];
fg   = {fg1};

normx1 = {@(x)1./(1+exp(-x)),@exp};
normx  = {normx1};

sim_sim_H(ii,simcat,simstr,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx,fg);
if ii(end)~=nsim, return; end
sim_wrap(nsim,simcat,simstr,modelnames,mnames,pnames,normx);
end

function run_b2(b,ii,nsim)

simcat      = 'sim_H1';

symname     = 'brl2';
simstr      = sprintf('%s[%+0.2f]',symname,b);

Nbar        = [20 20];

models      = {'model_rwneut1','model_rw2'};
tasks       = {'task_rwgo1','task_rwgo2'};
modelnames  = {'Biased RL','Dual-\alpha RL'};
mnames      = modelnames;
pnames      = { {'\alpha','\beta','\it b'},{'\alpha^+','\alpha^-','\beta'} };

tau1 = 4*[1 1 .25]';
tau2 = 4*[1 1 1]';
tau  = {tau1,tau2};

tx1  = [.1 1 b]';
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
simcats = repmat({'sim_H1','sim_H1','sim_H1','sim_H1'},2,1);
fitsimstrs(1,:) = {'fit_brl[+0.25].mat','fit_brl[+0.50].mat','fit_brl[+0.75].mat','fit_brl[+1.00].mat'}; 
fitsimstrs(2,:) = {'fit_brl2[+0.25].mat','fit_brl2[+0.50].mat','fit_brl2[+0.75].mat','fit_brl2[+1.00].mat'};                
xgroups = {'0.25','0.5','0.75','1'};
          
fsimfits  = fullfile(getdefaults('pipedir'),simcats,fitsimstrs);

[rate05,methods] = sim_stat_H1(fsimfits,0);

% close all;
plot_fig(xgroups,rate05,methods);
set(gcf,'name',mfilename);
end

function plot_fig(xgroups,power,methods)
[fs,fst,fsl,fsy,fsalpha,fsxt,fsA,fn,fnt,xsA,ysA,fpos0,siz0,colmap,alf,cmaphbi,bw,colmapsim]=fig_plot_properties; %#ok<*ASGLU>

fpos = [.3 fpos0(4)];
siz = fpos./fpos0(3:4).*siz0;

ysA = 1.05;

bwi = .7*bw;
%---------

nr = 2;
nc = 1;

figure;
set(gcf,'units','centimeters');
fsiz = get(gcf,'position');
fsiz(3:4) = siz;
fsiz(1:2) = fpos0(1:2).*fsiz(1:2);
set(gcf,'position',fsiz);
set(gcf,'units','normalized');

%--------
abc = 'ABCDEF';
for i=1:nr

    subplot(nr,nc, i);

    mx = power{i};
    mx(mx<.01) = .01;

    errorbarKxN(mx,0*mx,xgroups,'',colmap,0,bwi);
    alpha(gca,alf);
    set(gca,'fontsize',fs,'fontname',fn);
    ylim([0 1.2]);
    ytick = 0:.2:1;
    set(gca,'ytick',ytick);
    ylabel('Sensitivity','fontsize',fsy);
    if i==2
        xlabel('Effect size','fontsize',fsy);
    end

    title(sprintf('Performance of HBI t-test in scenario %d',i),'fontsize',fst,'fontname',fnt);
    
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname','Calibri');
    if i==1
        lg = legend(methods,'location','northwest','fontsize',fsl,'Orientation','vertical');
        set(lg,'linewidth',1);
    end
end
end
