function emp_2(i)
% Figure 14
% Analyzing PD and control empirical data using HBI

mode = '';
if nargin==0, mode = 'sumfig'; end
if nargin==1, mode = 'run'; end
getdefaults('addpath');

switch mode
    case 'run'
        run_perm(i);
    case 'sumfig'
        sumfig;
    otherwise
        error('Unknown mode: %s',mode);        
end

end

function p = permtest
nperm   = 1000;
ng      = 2;

gnames   = {'healthy','ON'};
permname = 'HON';

tempdir = getdefaults('tempdir');
fdir    = fullfile(tempdir,'emp_PD',permname);
pipedir = getdefaults('pipedir');
fsumfit = fullfile(pipedir,'emp_PD',permname,sprintf('fit_perm%d.mat',nperm)); makedir(fullfile(pipedir,'emp_PD',permname));
if ~exist(fsumfit,'file')
    for i=1:nperm
        for g=1:ng
            fname = fullfile(fdir,sprintf('hbmc_perm%04d_g%d.mat',i,g));
            fcbm = load(fname); cbm = fcbm.cbm;

            cbm = hbmc_exceedance(cbm);
            cbm = hbmc_errorbar(cbm);
            cbm = hbmc_output(cbm);        
            fit(i,g) = struct('hbi',cbm); %#ok<AGROW,NASGU>
        end
    end
    save(fsumfit,'fit');
end

sumfit = load(fsumfit); fit = sumfit.fit;

dnbar = nan(nperm,1);
for i=1:nperm
    n1 = fit(i,1).hbi.model_frequency/sum(fit(i,1).hbi.model_frequency)*[0 -1 1]';
    n2 = fit(i,2).hbi.model_frequency/sum(fit(i,2).hbi.model_frequency)*[0 -1 1]';
    dnbar(i,:) =  n2 - n1;
end

fnames{1} = fullfile(pipedir,'emp_PD',gnames{1},'hbi.mat');
fconfs{1} = fullfile(pipedir,'emp_PD',gnames{1},'config.mat');
fnames{2} = fullfile(pipedir,'emp_PD',gnames{2},'hbi.mat');
fconfs{2} = fullfile(pipedir,'emp_PD',gnames{2},'config.mat');

% stat_rt_2step(fnames{1})
[pxp,Nbar,mnames,mx,ex,pnames]= emp_stat(fnames,fconfs);
n1 = Nbar{1}*[0 1 -1]';
n2 = Nbar{2}*[0 1 -1]';

dn = n2 - n1;

p = mean(dn>dnbar);

pmin = 1/nperm;
p(p<pmin) = pmin;
end

function run_perm(ii)
if all(ii==0)
    emp_PD(0,'healthy');
    emp_PD(1,'healthy');
    emp_PD(2,'healthy');

    emp_PD(0,'ON');
    emp_PD(1,'ON');
    emp_PD(2,'ON');
end

nsim = 1000;
ii   = ii + (0:50:(nsim-1));

rng('shuffle');

emp_PD_perm(ii,{'healthy','ON'},'HON');
end

function sumfig
fsum  = fullfile(getdefaults('sumdir'),'emp_2','sum.mat');

if ~exist(fsum,'file')
    pipedir   = getdefaults('pipedir');
    fnames{1} = fullfile(pipedir,'emp_PD','ON','hbi.mat');
    fconfs{1} = fullfile(pipedir,'emp_PD','ON','config.mat');
    fnames{2} = fullfile(pipedir,'emp_PD','healthy','hbi.mat');
    fconfs{2} = fullfile(pipedir,'emp_PD','healthy','config.mat');

    % difference between the two groups using a permutation test
    p = permtest;
    [pxp,Nbar,mnames,mx,ex,pnames]= emp_stat(fnames,fconfs);
    mnames{1} = {'NL','RL','Dual-\alpha'};
    mnames{2} = {'NL','RL','Dual-\alpha'};
    pnames{2}{2} = '\alpha';
    save(fsum,'p','pxp','Nbar','mnames','mx','ex','pnames');
else
    sums    = load(fsum);
    p       = sums.p;
    pxp     = sums.pxp;
    Nbar    = sums.Nbar;
    mnames  = sums.mnames;
    mx      = sums.mx;
    ex      = sums.ex;
    pnames  = sums.pnames;
end


fprintf('difference between the two groups using permutation test:\n')
display(p);

% close all;
fig_plot(pxp,Nbar,mnames,mx,ex,pnames)
end

function fig_plot(pxp,Nbar,mnames,mx,ex,pnames)
pnames = sim_plot_adjust_pnames(pnames);

[fs,fst,fsl,fsy,fsalpha,fsxt,fsA,fn,fnt,xsA,ysA,fpos0,siz0,colmap,alf,cmaphbi,bw,colmapsim]=fig_plot_properties; %#ok<*ASGLU>

fpos = [.55 .5];
siz = fpos./fpos0(3:4).*siz0;

bw = bw*1.2;

nr = 2;
nc = 4;

figure;
set(gcf,'units','centimeters');
fsiz = get(gcf,'position');
fsiz(3:4) = siz;
fsiz(1:2) = fpos0(1:2).*fsiz(1:2);
set(gcf,'position',fsiz);
set(gcf,'units','normalized');

nf = length(pxp);

ki  = 1;
abc = 'ABCDEF';
dl  = [.02 .03];
rb  = [.9 .7];
%--
for f=1:nf
    pxpf    = pxp{f};
    Nbarf   = Nbar{f};
    mnamesf = mnames{f};
    pnamesf = pnames{f};
    mxf     = mx{f};
    exf     = ex{f};
    dlf     = dl(f);
    rbf     = rb(f);
    
    pxpf(pxpf<.01) = .01;
    Nbarf(Nbarf<.01) = .01;
    
    
    subplot(nr,nc,(f-1)*nc+1);    
    errorbarKxN(pxpf,pxpf*0,mnamesf,'',cmaphbi,0,bw);
    alpha(gca,alf);
    set(gca,'fontsize',fs,'fontname',fn);
    title('Model selection','fontsize',fst,'fontname',fnt);
    ylim([0 1.05]);
    ytick = 0:.2:1;
    set(gca,'ytick',ytick);
    ylabel('PXP','fontsize',fsy);
    ax = ancestor(gca, 'axes');
    xaxes = get(ax,'XAxis');
    set(xaxes,'Fontsize',fsxt);
    text(xsA,ysA,abc(f),'fontsize',fsA,'Unit','normalized','fontname','Calibri'); ki=ki+1;

    %--
    subplot(nr,nc,(f-1)*nc+2);    
    errorbarKxN(Nbarf,Nbarf*0,mnamesf,'',cmaphbi,0,bw);
    alpha(gca,alf);
    set(gca,'fontsize',fs,'fontname',fn);
    title('Model comparison','fontsize',fst,'fontname',fnt);
    ylim([0 1.05]);
    ytick = 0:.2:1;
    set(gca,'ytick',ytick);
    ylabel('Model frequency','fontsize',fsy);
    ax = ancestor(gca, 'axes');
    xaxes = get(ax,'XAxis');
    set(xaxes,'Fontsize',fsxt);

    %--
    hk = subplot(nr,nc,(f-1)*nc+(3:4)); 

    [~,kbest] = max(Nbarf);
    pos1 = get(hk,'position');
    xst = .5;
    yst = 1.07;
    text(xst,yst,sprintf('Parameters of %s',mnamesf{kbest}),'fontsize',fst,'fontname',fnt,...
        'Unit','normalized','fontweight','bold','Parent',hk,'HorizontalAlignment','center'); hold on;   
    set(hk,'visible','off');

    pos = pos1;
    np = size(mxf,2);

    x0 = pos1(1);    
    lp = 1/np* (pos1(3) -(np-1)*dlf);

    for i=1:np
        pos1 = pos;
        pos1(1) = x0+(i-1)*(lp+dlf);
        pos1(3) = (lp - dlf);

        axes('Position',pos1);    
        errorbarKxN(mxf(:,i),exf(:,i),pnamesf(i),'',cmaphbi,0,bw*rbf); hold on;
        alpha(gca,alf);
        set(gca,'fontsize',fs,'fontname',fn);

        ax = ancestor(gca, 'axes');
        xaxes = get(ax,'XAxis');
        set(xaxes,'Fontsize',fsalpha,'TickLabelInterpreter','latex');

        ytick = get(gca,'ytick');
        ytick(end) = [];
        ytick = ytick(1:2:end);
        set(gca,'ytick',ytick);
    end
end

end

