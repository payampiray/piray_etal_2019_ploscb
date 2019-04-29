function sim_sim_G2(ii,simcat,simstr,frandwalks,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx,simcat0,simstr0,cat)
rng('shuffle');
for i=ii
    compy_temp(simcat,simstr,i,simcat0,simstr0,models);
    simrun(simcat,simstr,i,frandwalks,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx,cat);
end
end

%-----------------------
function compy_temp(simcat,simstr,id,simcat0,simstr0,models)

simdir0  = fullfile(getdefaults('tempdir'),simcat0,simstr0,sprintf('sim%02d',id));
simdir   = fullfile(getdefaults('tempdir'),simcat,simstr,sprintf('sim%02d',id));

fsimsim0  = fullfile(simdir0,'sim.mat'); 
fdatasim0 = fullfile(simdir0,'data.mat');
fsimsim  = fullfile(simdir,'sim.mat'); 
fdatasim = fullfile(simdir,'data.mat');

mkdir(simdir);
if ~exist(fsimsim,'file')
copyfile(fsimsim0,fsimsim);
end
if ~exist(fdatasim,'file')
    fdata = load(fdatasim0);
    for m=1:3
    fdata.models{m} = str2func(models{m});
    end
    save(fdatasim,'-struct','fdata');
end
end

function simrun(simcat,simstr,id,frandwalks,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx,cat)

simdir  = fullfile(getdefaults('tempdir'),simcat,simstr,sprintf('sim%02d',id));
fprintf('sim directory is %s\n',simdir);

T = 201;

for k=1:length(models)
    models{k} = str2func(models{k});
    tasks{k} = str2func(tasks{k});
end

if ~exist(fullfile(simdir,'sim.mat'),'file')
simit(simdir,T,frandwalks,tasks,Nbar,mu,tau,models,modelnames,mnames,pnames,normx);
end
if exist(fullfile(simdir,'sim.mat'),'file')
run(simdir,cat);
end

end

function simit(simdir,T,frandwalks,tasks,Nbar,mu,tau,models,modelnames,mnames,pnames,normx)
N = sum(Nbar);
K = length(Nbar);

D    = nan(1,K);

for k=1:K
    D(k) = length(mu{k});
    if length(tau{k})~=D(k)
        error('!');
    end
end

beta = 1;

config  = struct('N',N,'K',K,'Nbar',Nbar,'D',D,'tau',{tau},'mu',{mu},...
    'beta',beta,'models',{models},'modelnames',{modelnames},'mnames',{mnames},'pnames',{pnames},'normx',{normx});

%---


z    = zeros(K,N);
data = cell(N,1);
h    = cell(1,K);
taskdata = cell(N,1);

init   = cell(K,1);
iz     = [];
for k=1:K
    init{k} = zeros(1,D(k));
    iz = [iz; k*ones(Nbar(k),1)];
end

for k=1:K
    t  = nan(D(k),N);
    t(:,iz==k)  = bsxfun(@times,((tau{k}).^-.5),randn(D(k),Nbar(k)));
    h{k}       = bsxfun(@plus,t, + mu{k} - nanmean(t,2));    
end
hk      = cell(1,K);
F    = nan(K,N);
for n=1:N      
    k = iz(n);
    theta   = h{k}(:,n);
    [rands]=randomwalk(T,frandwalks);    
    taskdata{n} = rands;
    [F(k,n),data{n}] = tasks{k}(rands,theta);
%     logpdf = @(x)models{k}(x,dataa);
%     x = slicesample(zeros(1,D(k)),2000,'logpdf',logpdf,'thin',5,'burnin',1000);      
end
for k=1:K
    nzk    = iz==k;
    z(k,nzk)  = 1;
end

sim    = struct('config',config,'taskdata',{taskdata},'z',z,'mu',{mu},'tau',{tau},'h',{h},'F',F); %#ok<NASGU>
% buildit(T,startpoint,endpoint,sdev,tasks,mu,tau,models,sim);

if exist(simdir,'dir')
%   rmdir(simdir,'s');
%     error('%s already exist',simdir);
end
makedir(simdir);

fconfig  = fullfile(simdir,'..','config.mat'); 
if ~exist(fconfig,'file')
    save(fconfig,'config');
end

fsimsim  = fullfile(simdir,'sim.mat'); 
fdatasim = fullfile(simdir,'data.mat');
if ~exist(fsimsim,'file')
    pause(ceil(5*rand));
    save(fsimsim,'sim');
end
if ~exist(fdatasim,'file')
    save(fdatasim,'data','models','init');
end


end

function rands = randomwalk(T,frandwalks)

rn = randperm(2,1);

walk     = load(frandwalks{rn});

rands.probs   = walk.probs';
rands.noise1  = rand(T,1);
rands.noise2  = rand(T,1);
rands.noiseR  = rand(T,1);
itr           = randperm(T, round(T/3) );
trans         = ones(T,1);
trans(itr)    = 2;
rands.trans   = trans;

end

%-----------------------
function [F] = model_hybrid(params,data)
% normx  = {@(x)(1./(1+exp(-x))),@(x)(1./(1+exp(-x))),@(x)(1./(1+exp(-x))),@exp,@(x)(1./(1+exp(-x))),@(x)x,@exp};
% pnames = {'\alpha_1','\alpha_2','\lambda','\beta_1','\it w','\it p','\beta_2'};

%-------------------------
fxu   = @(t)(1./(1+exp(-t)));
fxp   = @(t)exp(t);

alpha1  = fxu(params(1));
alpha2  = fxu(params(2));
lambda  = fxu(params(3));
tau     = fxp(params(4)); 
weight  = fxu(params(5));
phi     = params(6);
beta    = fxp(params(7));

%-------------------------
b     = [tau*weight tau*(1-weight) phi beta];

%%---
% unpack data
a1v = data.choice1;
a2v = data.choice2;
rv  = data.outcome;

missed = a1v==0 | a2v==0 | rv==0;
a1v(missed)=[];
a2v(missed)=[];
rv(missed)=[];
s2v = nan(length(a2v),1);
s2v(a2v==1 | a2v==2)=2;
s2v(a2v==3 | a2v==4)=3;
a2v(a2v==3 | a2v==4)=a2v(a2v==3 | a2v==4)-2;
rv (rv==2)=0;

T = length(rv);
n = zeros(3,2);
QTD = zeros(3,2);
QMB = zeros(1,2);
rep = zeros(1,2);

xQTD1 = nan(T,1);
xQTD2 = nan(T,1);
xQMB  = nan(T,1);
xrep  = nan(T,1);
for t=1:T
    s2 = s2v(t);
    a1 = a1v(t);
    a2 = a2v(t);
    
    nota1 = 3-a1;
    nota2 = 3-a2;
    xQTD1(t) = QTD(1,a1) - QTD(1,nota1);
    xQTD2(t) = QTD(s2,a2)- QTD(s2,nota2);
    xQMB (t) = QMB(1,a1) - QMB(1,nota1);
    xrep (t) = rep(1,a1) - rep(1,nota1);
    
    r  = rv(t);
    
    delta1 = 0 + QTD(s2,a2) - QTD(1,a1);
    QTD(1,a1) = QTD(1,a1) + alpha1*delta1;
    
    delta2 = r - QTD(s2,a2);
    QTD(s2,a2) = QTD(s2,a2) + alpha2*delta2;
    QTD(1,a1) = QTD(1,a1) + alpha1*lambda*delta2;
    
    % update transition probabilities
    n(s2,a1) = n(s2,a1)+1;    
    % transition to sB follwoing aA plus sC following aB
    n1 = n(2,1)+n(3,2);
    % or vice versa, to sC following aA plus sB following aB
    n2 = n(3,1)+n(2,2); 
    if n1>n2
        ptr(2,1) = 0.7;
        ptr(3,2) = 0.7;
        ptr(3,1) = 0.3;
        ptr(2,2) = 0.3;
    elseif n1<n2
        ptr(2,1) = 0.3;
        ptr(3,2) = 0.3;
        ptr(3,1) = 0.7;
        ptr(2,2) = 0.7;
    elseif n1==n2
        ptr(2,1) = 0.5;
        ptr(3,2) = 0.5;
        ptr(3,1) = 0.5;
        ptr(2,2) = 0.5;        
    end
    
    for aj=1:2
        QMB(aj) = ptr(2,aj)*max(QTD(2,:)) + ptr(3,aj)*max(QTD(3,:));
    end

    rep = rep*0;
    rep(a1) = 1;    
end

x      = [xQMB xQTD1 xrep];
X      = blkdiag(x,xQTD2);
%----

z = bsxfun(@times,X,b);
f = (1./(1+exp(-sum(z,2))));
F = sum(log(f+eps));
end

function [F] = model_mb(params,data)

fx = nan(1,7);
ip = [2 4 6 7];
fx(ip) = params;
fx(1) = -inf; % i.e alpha1=0
fx(3) = -inf; % lambda=0
fx(5) = inf;    % w=1

[F] = model_hybrid(fx,data);
end

function [F] = model_mf(params,data)

fx = nan(1,7);
ip = [1 2 3 4 6 7];
fx(ip) = params;
fx(5) = -inf;    % w=0

[F] = model_hybrid(fx,data);
end

%-----------------------
function [L,data]=task_hybrid(taskdata,params)

fxu   = @(t)(1./(1+exp(-t)));

alpha1  = fxu(params(1));
alpha2  = fxu(params(2));
lambda  = fxu(params(3));
tau     = exp(params(4)); 
omega   = fxu(params(5));
phi     = params(6);
beta    = exp(params(7));

tlin  = [tau omega phi beta]; %#ok<NASGU>
b     = [tau*omega tau*(1-omega) phi beta];

% unpack data
rprobs = taskdata.probs;
noise1 = taskdata.noise1;
noise2 = taskdata.noise2;
noiseR = taskdata.noiseR;
trans  = taskdata.trans;
T = length(noise1);
a1v = nan(T,1);
a2v = nan(T,1);
s2v = nan(T,1);
rv  = nan(T,1);    

T = length(a1v);
n = zeros(3,2);
QTD = zeros(3,2);
QMB = zeros(1,2);
rep = zeros(1,2);

p1 = nan(T,1);
p2 = nan(T,1);
for t=1:T
    
    
    
    xQTD1  = QTD(1,1) - QTD(1,2);
    xQMB   = QMB(1,1) - QMB(1,2);
    xrep   = rep(1,1) - rep(1,2);
    x      = [xQMB xQTD1 xrep];    
    p1(t)  = 1./(1+exp(-x*b(1:3)'));
    
%     pmb = 1./(1+exp(-x*[tau 0 phi]'));
%     pmf = 1./(1+exp(-x*[0 tau phi]'));
%     p1(t)  = omega*pmb + (1-omega)*pmf;
    
    %-------%-------
    % only for sim
        a1 = 2 - (noise1(t)<p1(t)); % i.e. if p>noise(t) -> a = 1; otherwise a = 2        
        s2 = a1+1;
        if trans(t)==2 % if rare
            s2 = 4-a1;
        end
        if s2==1
            disp('1');
        end
        
        a1v(t) = a1;
        s2v(t) = s2;
        f1(t,1)  = p1(t)*(a1==1) + (1-p1(t))*(a1==2);
    %-------%-------
    xQTD2  = QTD(s2,1) - QTD(s2,2);
    
    p2(t)  = (1./(1+exp(-beta*xQTD2)));
    
    %-------%-------
        a2     = 2 - (noise2(t)<p2(t)); % i.e. if p>noise(t) -> a = 1; otherwise a = 2        
        rv(t)   = 2 - (noiseR(t)<rprobs(t, (s2-2)*2 +a2  ));        
        a2v(t) = a2;
        f2(t,1)  = p2(t)*(a2==1) + (1-p2(t))*(a2==2);
    %-------%-------
    r  = (rv(t)==1)+0.0;
    
    delta1 = 0 + QTD(s2,a2) - QTD(1,a1);
    QTD(1,a1) = QTD(1,a1) + alpha1*delta1;
    
    delta2 = r - QTD(s2,a2);
    QTD(s2,a2) = QTD(s2,a2) + alpha2*delta2;
    QTD(1,a1) = QTD(1,a1) + alpha1*lambda*delta2;
    
    % update transition probabilities
    n(s2,a1) = n(s2,a1)+1;    
    % transition to sB follwoing aA plus sC following aB
    n1 = n(2,1)+n(3,2);
    % or vice versa, to sC following aA plus sB following aB
    n2 = n(3,1)+n(2,2); 
    if n1>n2
        ptr(2,1) = 0.7;
        ptr(3,2) = 0.7;
        ptr(3,1) = 0.3;
        ptr(2,2) = 0.3;
    elseif n1<n2
        ptr(2,1) = 0.3;
        ptr(3,2) = 0.3;
        ptr(3,1) = 0.7;
        ptr(2,2) = 0.7;
    elseif n1==n2
        ptr(2,1) = 0.5;
        ptr(3,2) = 0.5;
        ptr(3,1) = 0.5;
        ptr(2,2) = 0.5;        
    end
    
    for aj=1:2
        QMB(aj) = ptr(2,aj)*max(QTD(2,:)) + ptr(3,aj)*max(QTD(3,:));
    end

    rep = rep*0;
    rep(a1) = 1;    
end

loglik1 = log( f1);
loglik2 = log( f2);

L = sum(loglik1 + loglik2);

a2v(s2v==2) = a2v(s2v==2)+2;
data = struct('choice1',a1v,'transit',trans,'choice2',a2v,'outcome',rv);

end

function [L,data]=task_mb(taskdata,params)

fx = nan(1,7);
ip = [2 4 6 7];
fx(ip) = params;
fx(1) = -inf; % i.e alpha1=0
fx(3) = -inf; % lambda=0
fx(5) = inf;    % w=1

[L,data]=task_hybrid(taskdata,fx);
end

function [L,data]=task_mf(taskdata,params)

fx     = nan(1,7);
ip     = [1 2 3 4 6 7];
fx(ip) = params;
fx(5)  = -inf;    % w=0

[L,data]=task_hybrid(taskdata,fx);
end

%-----------------------
function run(simdir,kk)

[data,models,flap,init,v0,~,fsave,flog,fname] = sim_sim_dataload(simdir);

alg    = 'hierlap';
config1 = struct('verbose',1);

K = 3;
pconfig(1:K,1) = deal(cbm_config(length(init{1}),alg,[]));           
for k=1:K
    pconfig(k) = cbm_config(length(init{k}),alg,config1);
end
% config1.discard_bad = 1;

config2 = config1;
config2.save_data  = 0; 
config2.save_prog  = 0;
config1.numinit_up=100;
config1.tolgrad_liberal=1000;

nn = 1:length(data);
if kk<4
hbmc_lap(data,models(kk),flap(kk),init(kk),v0,config1,config2,nn);
end

if kk>3 && kk<7
hbmc_lap(data,models(kk-3),flap(kk-3),init(kk-3),v0,config1,config2);
end

if kk==7
if ~exist(fname,'file'), hbmc(data,models,flap,fname,pconfig,fsave,flog); end
end

if kk==8
if ~isempty(fsave), [fdir,fsave]=fileparts(fsave); fsave=fullfile(fdir,sprintf('%s0.mat',fsave)); end
if ~isempty(flog),  [fdir,flog]=fileparts(flog); flog=fullfile(fdir,sprintf('%s0.log',flog)); end
[fdir,fname]=fileparts(fname); fname=fullfile(fdir,sprintf('%s0.mat',fname));
if ~exist(fname,'file'), hbmc(data,models,flap,fname,pconfig,fsave,flog,1); end
end

end

