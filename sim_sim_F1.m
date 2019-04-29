function sim_sim_F1(ii,simcat,simstr,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx,taskmode)
if nargin<13, taskmode = 100; end
    
rng('shuffle');
for i=ii
    simrun(simcat,simstr,i,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx,taskmode);
end

end

function simrun(simcat,simstr,id,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx,taskmode)

simdir  = fullfile(getdefaults('tempdir'),simcat,simstr,sprintf('sim%02d',id));
fprintf('sim directory is %s\n',simdir);

[T,startpoint,endpoint,sdev]=task_profile(taskmode);

% buildit(T,startpoint,endpoint,sdev,tasks,mu,tau,models);

for k=1:length(models)
    models{k} = str2func(models{k});
    tasks{k} = str2func(tasks{k});
end

if ~exist(fullfile(simdir,'sim.mat'),'file')
simit(simdir,T,startpoint,endpoint,sdev,tasks,Nbar,mu,tau,models,modelnames,mnames,pnames,normx);
end
if exist(fullfile(simdir,'sim.mat'),'file')
run(simdir);
end
end

%-----------------------
function run(simdir)

[data,models,flap,init,v0,pconfig,fsave,flog,fname] = sim_sim_dataload(simdir);

config1.numinit_med = 30;
config1.numinit_up  = 30;
config1.discard_bad = 1;
config2.save_data  = 0; 
config2.save_prog  = 0;
nn = 1:length(data);
hbmc_lap(data,models,flap,init,v0,config1,config2,nn);
hbmc_lap(data,models,flap,init,v0,config1,config2);
if ~exist(fname,'file'), hbmc(data,models,flap,fname,pconfig,fsave,flog); end
if ~isempty(fsave), [fdir,fsave]=fileparts(fsave); fsave=fullfile(fdir,sprintf('%s0.mat',fsave)); end
if ~isempty(flog),  [fdir,flog]=fileparts(flog); flog=fullfile(fdir,sprintf('%s0.log',flog)); end
[fdir,fname]=fileparts(fname); fname=fullfile(fdir,sprintf('%s0.mat',fname));
if ~exist(fname,'file'), hbmc(data,models,flap,fname,pconfig,fsave,flog,1); end
end

%--------------------------
function simit(simdir,T,startpoint,endpoint,sdev,tasks,Nbar,mu,tau,models,modelnames,mnames,pnames,normx)
N = sum(Nbar);
K = length(Nbar);

if K~=length(models)
    error('!');
end

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
    [r]=randomwalk(T,startpoint,endpoint,sdev);    
    [R,taskdata{n}] = tasks{k}(r,theta);
    dataa = taskdata{n};
    dataa.outcome = R;
    [F(k,n),data{n}] = models{k}(theta,dataa);    
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
    save(fsimsim,'sim');
end

% This is the only change with respect to simit in sim_sim
models{2} = @model_rw0;
init{2}   = 0;
models(3:end) = []; %#ok<NASGU>
init(3:end) = []; %#ok<NASGU>
if ~exist(fdatasim,'file')
    save(fdatasim,'data','models','init');
end


end
%--------------------------
function [loglik,data,p] = model_rw0(params,dat)
% fxu     = @(x)1./(1+exp(-x));
alpha   = 1;
beta    = exp(params(1));

simmode = 0;
if ~isfield(dat,'actions'), simmode = 1; end

if simmode, noise = dat.noise; end

% unpack data
R         = dat.outcome;
T         = size(R,1);

epsilon = 0;
if simmode==0, A = dat.actions; epsilon = 0; end

q       = zeros(1,2);
f       = nan(T,1);
actions = nan(T,1);
outcome = nan(T,1);

p = nan(T,1);
n = [0 0];
for t=1:T    
    p1     = 1./(1+exp(-beta*(q(1)-q(2))));
    p(t)   = epsilon/2 + (1-epsilon)*p1;
    
    if simmode     
        a    = 2 - (noise(t)<p(t)); % i.e. if p>noise(t) -> a = 1; otherwise a = 2
        n(a) = n(a)+1;
        o    = R(n(a),a);
    else
        a    = A(t);
        o    = R(t);
    end    
    f(t)     = p(t).*(a==1) + (1-p(t)).*(a==2);
    
    delta    = o - q(a)';
    q(a)     = q(a) + (alpha.*delta)';
    
    actions(t) = a;
    outcome(t) = o;
    
end

if simmode
    data = struct('actions',actions,'outcome',outcome,'simmed',dat);
end

loglik = sum(log(f+eps));
end

function [loglik,data,p] = model_rw1(params,dat)
fxu     = @(x)1./(1+exp(-x));
alpha   = fxu(params(1));
beta    = exp(params(2));

simmode = 0;
if ~isfield(dat,'actions'), simmode = 1; end

if simmode, noise = dat.noise; end

% unpack data
R         = dat.outcome;
T         = size(R,1);

epsilon = 0;
if simmode==0, A = dat.actions; epsilon = 0; end

q       = zeros(1,2);
f       = nan(T,1);
actions = nan(T,1);
outcome = nan(T,1);

p = nan(T,1);
n = [0 0];
for t=1:T    
    p1     = 1./(1+exp(-beta*(q(1)-q(2))));
    p(t)   = epsilon/2 + (1-epsilon)*p1;
    
    if simmode     
        a    = 2 - (noise(t)<p(t)); % i.e. if p>noise(t) -> a = 1; otherwise a = 2
        n(a) = n(a)+1;
        o    = R(n(a),a);
    else
        a    = A(t);
        o    = R(t);
    end    
    f(t)     = p(t).*(a==1) + (1-p(t)).*(a==2);
    
    delta    = o - q(a)';
    q(a)     = q(a) + (alpha.*delta)';
    
    actions(t) = a;
    outcome(t) = o;
    
end

if simmode
    data = struct('actions',actions,'outcome',outcome,'simmed',dat);
end

loglik = sum(log(f+eps));
end
%--------------------------
function [R,data]=task_rw0(r,params) %#ok<INUSD>
% fxu     = @(x)1./(1+exp(-x));
alpha1  = 0;

T = length(r);

x = [0 0];
X = nan(T,2);
for t=1:T
    delta1 = (r(t,:)-x);
%     delta2 = (r(t)-y);
%     alpha  = (delta2>=0)*alpha1 + (delta2<0)*alpha2;
    x = x + alpha1*delta1;
    X(t,:) = x;    
%     y = y + alpha*delta2;
end

R = X;
R = 2*(X>0)-1;

% for i=1:10
%     c(i) = corr(X((i+1):end),X(1:end-i));
% end

data.outcome = R;
data.noise   = rand(T,1);
data.r       = r;
data.X       = [];
end

function [R,data]=task_rw1(r,params)
fxu     = @(x)1./(1+exp(-x));
alpha1  = fxu(params(1));

T = length(r);

x = [0 0];
X = nan(T,2);
for t=1:T
    delta1 = (r(t,:)-x);
%     delta2 = (r(t)-y);
%     alpha  = (delta2>=0)*alpha1 + (delta2<0)*alpha2;
    x = x + alpha1*delta1;
    X(t,:) = x;    
%     y = y + alpha*delta2;
end

R = X;
R = 2*(X>0)-1;

% for i=1:10
%     c(i) = corr(X((i+1):end),X(1:end-i));
% end

data.outcome = R;
data.noise   = rand(T,1);
data.r       = r;
data.X       = [];
end
%--------------------------
function [T,startpoint,endpoint,sdev]=task_profile(mode)
if nargin<1, mode =3; end

switch mode
    case 100
        T=100; 
    case 50
        T=50; 
    case 200
        T=200; 
    otherwise
        error('Unknown mode!');
end

startpoint = -1;
endpoint   = +1;
sdev = 1;
end

function [r,walk]=randomwalk(N,startpoint,endpoint,s)
steps = randn(N,2)*s;
% ensure the sum is exactly zero
steps = bsxfun(@minus,steps , mean(steps));
% add in a bias to each step.
ds = (endpoint - startpoint);
steps = bsxfun(@plus,steps , [ds -ds]/N);
walk = cumsum(steps);

r = steps;
end