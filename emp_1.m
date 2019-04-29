function emp_1(i)
% Figure 13
% Analyzing 2-step empirical data using HBI

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
emp_2step(ii)
end

function p = stat_rt_2step(fname)
pipedir  = getdefaults('pipedir');
frt = fullfile(pipedir,'emp_2step','sum_RT.mat');

if ~exist(frt,'file')
    RT      = load(fullfile(pipedir,'emp_2step','2step_RT.mat')); RT = RT.RT;

    N = length(RT);
    rt = nan(N,1);
    for n=1:N
        rt(n,:) = nanmedian(RT{n}(:,1) )/1000;
    end
    save(frt,'rt');
end
frt = load(frt); rt = frt.rt;



cbm  = load(fname); cbm = cbm.cbm;
r   = cbm.r;
[~,rm] = max(r,[],1);
if any(rm>2)
    error('!');
end

%--------------------------------------------------------------------
mrt = nan(1,2);
srt = nan(1,2);
for m=1:2
    mrt(m) = mean(rt(rm==m));
    srt(m) = serr(rt(rm==m));
end

ii = [1 2];
p  = ranksum(rt(rm==ii(1)),rt(rm==ii(2)));

end

function sumfig

fsum  = fullfile(getdefaults('sumdir'),'emp_1','sum.mat');

if ~exist(fsum,'file')
    pipedir = getdefaults('pipedir');
    fnames{1} = fullfile(pipedir,'emp_2step','hbi.mat');
    fconfs{1} = fullfile(pipedir,'emp_2step','config.mat');
    p = stat_rt_2step(fnames{1});
    [pxp,Nbar,mnames,mx,ex,pnames]= emp_stat(fnames,fconfs);
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

fprintf('difference in RT between the two subgroups:\n')
display(p);


% close all;
fig_plot(pxp,Nbar,mnames,mx,ex,pnames)
end

function fig_plot(pxp,Nbar,mnames,mx,ex,pnames)
pnames = sim_plot_adjust_pnames(pnames);

[fs,fst,fsl,fsy,fsalpha,fsxt,fsA,fn,fnt,xsA,ysA,fpos0,siz0,colmap,alf,cmaphbi,bw,colmapsim]=fig_plot_properties; %#ok<*ASGLU>

fpos = [.6 .27];
siz = fpos./fpos0(3:4).*siz0;
% siz = siz0;

bw = bw*1.2;
posm = 220; % in pixel
nr = 1;
nc = 4;

figure;
set(gcf,'units','centimeters');
fsiz = get(gcf,'position');
fsiz(3:4) = siz;
fsiz(1:2) = fpos0(1:2).*fsiz(1:2);
set(gcf,'position',fsiz);
set(gcf,'units','normalized');

ki  = 1;
abc = 'ABCDEF';
dl  = .015;
rb  = 2.5;

%--
f = 1;
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
set(gca,'units','pixels');
posp = get(gca,'position');
posp(4) = posm;
posp(2) = 1.2*posp(2);
set(gca,'position',posp);    


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
text(xsA,ysA,abc(1),'fontsize',fsA,'Unit','normalized','fontname','Calibri'); ki=ki+1;

%--
subplot(nr,nc,(f-1)*nc+2);    
set(gca,'units','pixels');
posp = get(gca,'position');
posp(4) = posm;
posp(2) = 1.2*posp(2);
set(gca,'position',posp);     

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
text(xsA,ysA,abc(2),'fontsize',fsA,'Unit','normalized','fontname','Calibri');

%--
hk = subplot(nr,nc,(f-1)*nc+(3:4)); 
set(gca,'units','pixels');
posp = get(gca,'position');
posp(4) = posm;
posp(2) = 1.2*posp(2);
set(gca,'position',posp);   
set(gca,'units','normalized');

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

    if i==1
        text(xsA-1,ysA,'C','fontsize',fsA,'Unit','normalized','fontname','Calibri','HorizontalAlignment','left');
    end

    ytick = get(gca,'ytick');
    ytick(end) = [];
    ytick = ytick(1:2:end);
    set(gca,'ytick',ytick);


end

end

