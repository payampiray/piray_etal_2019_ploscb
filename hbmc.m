function math = hbmc(data,models,flap,fname,pconfig,fsave,flog,isnull)
if nargin<8, isnull = 0; end
load_hbmc = [];

%==========================================================================
a0 = 0;
b = 1; v = 0.5; s = 0.01;
prior = struct('b',b,'v',v,'s',s);

%==========================================================================
N = size(data,1);
K = length(models);

modelnames = cell(size(models));
for k=1:K
    modelnames{k} = func2str(models{k});
end

fid = 1;
if ~usejava('desktop') && ~isempty(flog), fid = fopen(flog,'w'); end

input = struct('data',[],'models',{models},'flap',{flap},'config',pconfig);
%--------------------------------------------------------------------------
% log and report
clc;
verbose = 10;

if verbose, fprintf(fid, '%-40s%30s\n',mfilename,datestr(now)); end
if verbose, fprintf(fid, '%-70s\n',repmat('=',1,70)); end

%--------------------------------------------------------------------------
[cm,thetabar,Sdiag,pmutau,pm,bound,prog] = hbmc_init(input.flap,a0,b,v,s,isnull);
Nbar = N*ones(K,1); r = ones(K,N);
dprog.dalphaprog = [];
optim.mutau = struct('optim_v',0,'optim_b',0);
optim.cm.config = pconfig;
iter0 = 1;

if ~isempty(load_hbmc)
    math0 = load(load_hbmc); math0 = math0.math;
    fsave = load_hbmc;
    
    iter0 = length(math0)+1;
    iter  = length(math0);
    r     = math0(iter).r;    
    Nbar  = math0(iter).Nbar;    
    thetabar = math0(iter).thetabar;
    Sdiag    = math0(iter).Sdiag;    
    qmutau   = math0(iter).qmutau;
    pmutau   = math0(iter).pmutau;
    qm       = math0(iter).qm;
    pm       = math0(iter).pm;
    
    optim = math0(iter).optim;
    bound = math0(iter).bound;
    prog  = math0(iter).prog;
    dprog = math0(iter).dprog;    
end

maxiter = 50;
for iter = iter0:maxiter    
    fprintf(fid, '%s\n',repmat('-',1,70));        
    fprintf(fid, 'Iteration %02d\n',iter);
    
    if iter>1
        fprintf(fid, '-- Quadratic approximation\n');
        [cm] = hbmc_cm(models,data,pconfig,qmutau,cm);

        fprintf(fid, '-- Computing q(H,Z)\n');
        [r,Nbar,thetabar,Sdiag,bound.qHZ] = hbmc_qHZ(qmutau,qm,cm);
        ss = '';
        for k=1:K, ss = sprintf('%s%s: %2.1f|',ss,modelnames{k},Nbar(k)); end
        fprintf(fid,'\n-----%s\n',ss);        
        
        [bound,dL] = hbmc_bound(bound,'qHZ');
        fprintf(fid,'%-40s%30s\n',' ',sprintf('dL: %7.5f',dL));
    end
    
    fprintf(fid, '-- Computing q(mu,tau)\n');
    [qmutau,pmutau,bound.qmutau,optim.mutau] = hbmc_qmutau(pmutau,Nbar,thetabar,Sdiag,optim.mutau);
    [bound,dL] = hbmc_bound(bound,'qmutau');
    fprintf(fid,'%-40s%30s\n',' ',sprintf('dL: %7.5f',dL));
    
    fprintf(fid, '-- Computing q(m)\n');
    [qm,bound.qm] = hbmc_qm(pm,Nbar);
    [bound,dL] = hbmc_bound(bound,'qm');
    fprintf(fid,'%-40s%30s\n',' ',sprintf('dL: %7.5f',dL));
    
    
    [dprog,prog] = cal_prog(prog,bound.bound.L,r,qm.alpha,thetabar,Sdiag,dprog.dalphaprog);
    math(iter)   = struct('cm',cm,...
                          'r',r,...
                          'Nbar',Nbar,'thetabar',{thetabar},'Sdiag',{Sdiag},...
                          'pm',pm,'pmutau',pmutau,'qm',qm,'qmutau',qmutau,...                          
                          'bound',bound,'prog',prog,'dprog',dprog,...
                          'optim',optim,'input',input,'prior',prior); %#ok<AGROW>
    if ~isempty(fsave), save(fsave,'math'); end
    
    fprintf(fid,'\nIteration %02d summary:\n',iter);
    fprintf(fid,'%-40s%30s\n',' ',sprintf('dL: %7.2f',dprog.dL));
    fprintf(fid,'%-40s%30s\n',' ',sprintf('dN: %7.2f',dprog.dalpha));
    fprintf(fid,'%-40s%30s\n',' ',sprintf('dx: %7.2f',dprog.dx));
    if terminate(dprog,0)
        fprintf(fid,'%-40s%30s\n',' ','Converged!');
        break;
    end
end

if ~isempty(fname)
save(fname,'math');
end
if fid~=1, fclose(fid); end
end

function done = terminate(dprog,mode)
if nargin<2, mode = 1; end

toldalpha = 0.05;
toldx     = 0.01;
dalpha     = dprog.dalpha;
dx         = dprog.dx;

switch mode
    case 1
        done    = dalpha<toldalpha; %(dalphaprog < toldalpha);
    case 0
        done    = dx<toldx;
end
end

function [dprog,prog] = cal_prog(prog,L,r,alpha,thetabar,Sdiag,dalphaprog)
i         = length(prog);
if ~isempty(prog)
    Lpre      = prog(i).L;
    rpre      = prog(i).r;
    alphapre  = prog(i).alpha;
    xpre      = prog(i).x;
else
    xpre      = nan;
    Lpre      = nan;
    rpre      = nan(size(r));
    alphapre  = nan(size(alpha));
    prog      = struct('L',L,'r',r,'alpha',alpha,'x',nan);
end

K = length(thetabar);
x = [];
for k=1:K
    x = [x; thetabar{k}./sqrt(Sdiag{k})]; %#ok<AGROW>
end

dx            = sqrt(mean((x-xpre).^2));
dL            = L-Lpre;
dr            = mean(mean(abs(r-rpre)));

[~,ibest]     = max(alpha);
dalpha        = abs(alpha(ibest)-alphapre(ibest));
dalphaprog    = [dalphaprog dalpha];
dprog         = struct('dL',dL,'dr',dr,'dalpha',dalpha,'dalphaprog',dalphaprog,'dx',dx);
prog(i+1)     = struct('L',L,'r',r,'alpha',alpha,'x',x);
end
