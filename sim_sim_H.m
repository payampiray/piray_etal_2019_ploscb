function sim_sim_H(ii,simcat,simstr,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx,free_group,anormal)
if nargin<14, anormal = []; end
    
rng('shuffle');
for i=ii
    simrun(simcat,simstr,i,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx,free_group,anormal);
end

end

function simrun(simcat,simstr,id,Nbar,models,tasks,modelnames,mnames,pnames,mu,tau,normx,free_group,anormal)
do_anormal = 0;
if ~isempty(anormal), do_anormal=1; end

simdir  = fullfile(getdefaults('tempdir'),simcat,simstr,sprintf('sim%02d',id));
fprintf('sim directory is %s\n',simdir);

[T,startpoint,endpoint,sdev,nS]=task_profile;

models0 = models;
if exist(fullfile(simdir,'sim.mat'),'file')
    fsimsim  = fullfile(simdir,'sim.mat'); 
    fdatasim = fullfile(simdir,'data.mat');
    
    sim = load(fsimsim); sim=sim.sim;
    sim.models = models0;
    sim.config.models = models0;
    save(fsimsim,'sim');
    
    fdata = load(fdatasim);
    fdata.models = models0;
    save(fdatasim,'-struct','fdata');    
end

if ~exist(fullfile(simdir,'sim.mat'),'file')
    if ~do_anormal        
        simit(simdir,T,startpoint,endpoint,sdev,nS,tasks,Nbar,mu,tau,models,modelnames,mnames,pnames,normx);
    else        
        simanit(simdir,T,startpoint,endpoint,sdev,nS,tasks,Nbar,mu,tau,models,modelnames,mnames,pnames,normx,anormal);
    end
end
if exist(fullfile(simdir,'sim.mat'),'file')
run0(simdir,free_group);
end
end

%--------------------------
function [loglik,data,p,dq] = model_rwneut1(params,dat)
fxu    = @(x)1./(1+exp(-x));
alpha  = fxu(params(1));
beta   = exp(params(2));
c      = params(3);
b      = 0;

gobias = [1 -1 0];

simmode = 0;
if ~isfield(dat,'actions'), simmode = 1; end

if simmode, noise = dat.noise; end

% unpack data
R         = dat.outcome;
S         = dat.stimulus;
T         = size(R,1);

nS = max(S);

epsilon = 0;
if simmode==0, A = dat.actions; epsilon = 0; end

q       = zeros(nS,2); % stimulus x action
f       = nan(T,1);
actions = nan(T,1);
outcome = nan(T,1);
dq      = nan(T,1);
p = nan(T,1);

n = [0 0];
for t=1:T
    s      = S(t);
    
    p1     = 1./(1+exp(-beta*(q(s,1)-q(s,2))-(c+b)*gobias(s) ));
    p(t)   = epsilon/2 + (1-epsilon)*p1;
    
    dq(t)      = q(s,1)-q(s,2);
    if simmode     
        a    = 2 - (noise(t)<p(t)); % i.e. if p>noise(t) -> a = 1; otherwise a = 2
        n(a) = n(a)+1;
        o    = R(n(a),a);
    else
        a    = A(t);
        o    = R(t);
    end    
    f(t)     = p(t).*(a==1) + (1-p(t)).*(a==2);
    
    delta    = o - q(s,a)';
    q(s,a)   = q(s,a) + (alpha.*delta)';
    
    actions(t) = a;
    outcome(t) = o;
end

if simmode
    data = struct('stimulus',S,'actions',actions,'outcome',outcome,'simmed',dat);
end

loglik = sum(log(f+eps));
end

function [loglik,data,p,dq] = model_rw2(params,dat)
fxu    = @(x)1./(1+exp(-x));
alpha1 = fxu(params(1));
alpha2 = fxu(params(2));
beta   = exp(params(3));
c      = 0; %params(4);

gobias = [1 -1 0];

simmode = 0;
if ~isfield(dat,'actions'), simmode = 1; end

if simmode, noise = dat.noise; end

% unpack data
R         = dat.outcome;
S         = dat.stimulus; % 2 stimuli
T         = size(R,1);

nS = max(S);

epsilon = 0;
if simmode==0, A = dat.actions; epsilon = 0; end

q       = zeros(nS,2); % stimulus x action
f       = nan(T,1);
actions = nan(T,1);
outcome = nan(T,1);

p  = nan(T,1);
dq = nan(T,1);

n = [0 0];
for t=1:T
    s      = S(t);
    
    p1     = 1./(1+exp(-beta*(q(s,1)-q(s,2))-c*gobias(s) ));
    p(t)   = epsilon/2 + (1-epsilon)*p1;
    
    dq(t) = q(s,1)-q(s,2);
    if simmode     
        a    = 2 - (noise(t)<p(t)); % i.e. if p>noise(t) -> a = 1; otherwise a = 2
        n(a) = n(a)+1;
        o    = R(n(a),a);
    else
        a    = A(t);
        o    = R(t);
    end    
    f(t)     = p(t).*(a==1) + (1-p(t)).*(a==2);
    
    delta    = o - q(s,a)';
    alpha  = (delta>=0)*alpha1 + (delta<0)*alpha2;    
    q(s,a)   = q(s,a) + (alpha.*delta)';
    
    actions(t) = a;
    outcome(t) = o;
    
end

if simmode
    data = struct('stimulus',S,'actions',actions,'outcome',outcome,'simmed',dat);
end

loglik = sum(log(f+eps));
end

%--------------------------
function [R,data]=task_rwgo1(r,params, model)
S  = r.S;
r  = r.r;
nS = max(S);

fxu     = @(x)1./(1+exp(-x));
alpha1  = fxu(params(1));

T = length(r);

x = zeros(nS,2);
X = nan(T,2);
for t=1:T
    s      = S(t);
    delta1 = (r(t,:)-x(s,:));
    x(s,:) = x(s,:) + alpha1*delta1;
    X(t,:) = x(s,:);
end

R = 2*(X>0)-1;

% make sure there is no way to get reward by having a biased go approach
dR = R*[.5 -.5]';
for s=1:2
    Rd    = R(S==s,:)*[.5 -.5]';
    nbias = sum(Rd);
    ii = find(S==s & dR==sign(nbias));
    tts = randperm(length(ii));
    if mod(abs(nbias),2)==1
        R(ii(tts(1)),:) = [0 0];
        tts(1) = [];
    end
    
    tts = ii(tts(1:floor(abs(nbias)/2)));
    rr = R(tts,:);
    R(tts,1) = rr(:,2);
    R(tts,2) = rr(:,1);
    
    ncorrbiased(s) = sum(R(S==s,:)*[.5 -.5]'); %#ok<AGROW,NASGU> % should be zero
end
    
data.stimulus= S;
data.outcome = R;
data.r       = r;

[data,err,i] = noise_analysis(S,model,params,data);

end

function [R,data]=task_rwgo2(r,params,model)
S  = r.S;
r  = r.r;
nS = max(S);

fxu     = @(x)1./(1+exp(-x));
alpha1  = fxu(params(1));
alpha2  = fxu(params(2));

T = length(r);

x = zeros(nS,2);
X = nan(T,2);
for t=1:T
    s      = S(t);
    delta  = (r(t,:)-x(s,:));
    alpha  = (delta>=0)*alpha1 + (delta<0)*alpha2;
    
    x(s,:) = x(s,:) + alpha.*delta;
    X(t,:) = x(s,:);
end

R = 2*(X>0)-1;

% make sure there is no way to get reward by having a biased go approach
dR = R*[.5 -.5]';
for s=1:nS
    Rd    = R(S==s,:)*[.5 -.5]';
    nbias = sum(Rd);
    ii = find(S==s & dR==sign(nbias));
    tts = randperm(length(ii));
    if mod(abs(nbias),2)==1
        R(ii(tts(1)),:) = [0 0];
        tts(1) = [];
    end
    
    tts = ii(tts(1:floor(abs(nbias)/2)));
    rr = R(tts,:);
    R(tts,1) = rr(:,2);
    R(tts,2) = rr(:,1);
    
    ncorrbiased(s) = sum(R(S==s,:)*[.5 -.5]'); %#ok<AGROW,NASGU> % should be zero
end

data.stimulus= S;
data.outcome = R;
data.noise   = rand(T,1);
data.r       = r;

[data,err,i] = noise_analysis(S,model,params,data);
end

%-----------------------
function run0(simdir,free_group)
m0 = 0;
[data,models,flap,init,v0,pconfig,fsave,flog,fname] = sim_sim_dataload(simdir);

for k=1:length(models)
    models{k} = str2func(models{k});
end

config1.numinit_med = 30;
config1.numinit_up  = 30;
config1.discard_bad = 1;
config2.save_data  = 0;
config2.save_prog  = 0;
nn = 1:length(data);
hbmc_lap(data,models,flap,init,v0,config1,config2,nn);
hbmc_lap(data,models,flap,init,v0,config1,config2);
% hierlap_fixedgroup(data,models,flap,config2,free_group);
hierlap_fixedmean(data,models,flap,config2,free_group);

if ~exist(fname,'file'), hbmc(data,models,flap,fname,pconfig,fsave,flog); end
if ~isempty(fsave), [fdir,fsave]=fileparts(fsave); fsave=fullfile(fdir,sprintf('%s0.mat',fsave)); end
if ~isempty(flog),  [fdir,flog]=fileparts(flog); flog=fullfile(fdir,sprintf('%s0.log',flog)); end
[fdir,fname]=fileparts(fname); fname=fullfile(fdir,sprintf('%s0.mat',fname));
if ~exist(fname,'file'), hbmc(data,models,flap,fname,pconfig,fsave,flog,1); end
end

%--------------------------
function simit(simdir,T,startpoint,endpoint,sdev,nS,tasks,Nbar,mu,tau,funcmodels,modelnames,mnames,pnames,normx)
N = sum(Nbar);
K = length(Nbar);

D    = nan(1,K);

models = funcmodels;
for k=1:K
    funcmodels{k} = str2func(funcmodels{k});
    tasks{k} = str2func(tasks{k});    
    
    D(k) = length(mu{k});
    if length(tau{k})~=D(k)
        error('!');
    end
end

beta = 1;

config  = struct('N',N,'K',K,'Nbar',Nbar,'D',D,'tau',{tau},'mu',{mu},'beta',beta,'models',{models},...
                 'modelnames',{modelnames},'mnames',{mnames},'pnames',{pnames},'normx',{normx});

%---


z    = zeros(K,N);
data = cell(N,1);
h    = cell(1,K);
taskdata = cell(N,1);

init   = cell(K,1);
iz     = [];
for k=1:K
    init{k} = zeros(1,D(k));
    iz = [iz; k*ones(Nbar(k),1)]; %#ok<AGROW>
end

for k=1:K
    t  = nan(D(k),N);
    t(:,iz==k)  = bsxfun(@times,((tau{k}).^-.5),randn(D(k),Nbar(k)));

    % -- normalize
%     t = bsxfun(@minus, t, nanmean(t,2));
    % --
    
    h{k}        = bsxfun(@plus,t, + mu{k} );
end
F    = nan(K,N);
for n=1:N      
    k = iz(n);
    theta   = h{k}(:,n);
    [r]=randomwalk(T,startpoint,endpoint,sdev,nS);    
    [R,taskdata{n}] = tasks{k}(r,theta,funcmodels{k});
%     taskdata{n}.X.noisesim
    dataa = taskdata{n};
    dataa.outcome = R;
    [F(k,n),data{n}] = funcmodels{k}(theta,dataa);    
%     logpdf = @(x)models{k}(x,dataa);
%     x = slicesample(zeros(1,D(k)),2000,'logpdf',logpdf,'thin',5,'burnin',1000);      
end
for k=1:K
    nzk    = iz==k;
    z(k,nzk)  = 1;
end

sim    = struct('config',config,'taskdata',{taskdata},'z',z,'mu',{mu},'tau',{tau},'h',{h},'F',F); %#ok<NASGU>

if exist(simdir,'dir')
%   rmdir(simdir,'s');
%     error('%s already exist',simdir);
end
makedir(simdir);

fsimsim  = fullfile(simdir,'sim.mat'); 
fdatasim = fullfile(simdir,'data.mat');
if ~exist(fsimsim,'file')
    save(fsimsim,'sim');
end
if ~exist(fdatasim,'file')
    save(fdatasim,'data','models','init');
end


end

function simanit(simdir,T,startpoint,endpoint,sdev,nS,tasks,Nbar,mu,tau,funcmodels,modelnames,mnames,pnames,normx,anormal)
N = sum(Nbar);
K = length(Nbar);

D    = nan(1,K);

models = funcmodels;
for k=1:K
    funcmodels{k} = str2func(funcmodels{k});
    tasks{k} = str2func(tasks{k});    
    
    D(k) = length(mu{k});
    if length(tau{k})~=D(k)
        error('!');
    end
end

beta = 1;

config  = struct('N',N,'K',K,'Nbar',Nbar,'D',D,'tau',{tau},'mu',{mu},'beta',beta,'models',{models},...
                 'modelnames',{modelnames},'mnames',{mnames},'pnames',{pnames},'normx',{normx},'anormal',{anormal});

%---


z    = zeros(K,N);
data = cell(N,1);
h    = cell(1,K);
taskdata = cell(N,1);

init   = cell(K,1);
iz     = [];
for k=1:K
    init{k} = zeros(1,D(k));
    iz = [iz; k*ones(Nbar(k),1)]; %#ok<AGROW>
end

for k=1:K
    t  = nan(D(k),N);    
    if isempty(anormal{k})
        t(:,iz==k)  = bsxfun(@times,((tau{k}).^-.5),randn(D(k),Nbar(k)));    
        t           = bsxfun(@plus,t, + mu{k} );
    else
        for j=1:D(k)
            frnd_kj = anormal{k}{j};
            if isempty(frnd_kj)
                t(j,iz==k)  = bsxfun(@times,((tau{k}(j)).^-.5),randn(1,Nbar(k)));
                t(j,:)      = bsxfun(@plus,t(j,:), + mu{k}(j) );                
            else
                t(j,iz==k)  = frnd_kj(Nbar(k));                
            end
        end
    end
    
    h{k} = t;
    
end
F    = nan(K,N);
for n=1:N      
    k = iz(n);
    theta   = h{k}(:,n);
    [r]=randomwalk(T,startpoint,endpoint,sdev,nS);    
    [R,taskdata{n}] = tasks{k}(r,theta,funcmodels{k});
%     taskdata{n}.X.noisesim
    dataa = taskdata{n};
    dataa.outcome = R;
    [F(k,n),data{n}] = funcmodels{k}(theta,dataa);    
%     logpdf = @(x)models{k}(x,dataa);
%     x = slicesample(zeros(1,D(k)),2000,'logpdf',logpdf,'thin',5,'burnin',1000);      
end
for k=1:K
    nzk    = iz==k;
    z(k,nzk)  = 1;
end

sim    = struct('config',config,'taskdata',{taskdata},'z',z,'mu',{mu},'tau',{tau},'h',{h},'F',F); %#ok<NASGU>

if exist(simdir,'dir')
%   rmdir(simdir,'s');
%     error('%s already exist',simdir);
end
makedir(simdir);

fsimsim  = fullfile(simdir,'sim.mat'); 
fdatasim = fullfile(simdir,'data.mat');
if ~exist(fsimsim,'file')
    save(fsimsim,'sim');
end
if ~exist(fdatasim,'file')
    save(fdatasim,'data','models','init');
end


end

%--------------------------
function [data,err,i] = noise_analysis(S,model,params,data)
T = length(S);
beta = exp(params(2));
c = (params(3));


err = 1; i = 0; tolerr = 0.01; maxnoisesim = 5000;
while (err>tolerr) && i<maxnoisesim
    i=i+1;
    noise{i} = rand(T,1); %#ok<AGROW>
    data.noise   = noise{i};
    [~,datx,~,dv] = model(params,data);
    x=[dv -2*S+3];

    ex(:,i) = mean(prod(x,2)); %#ok<AGROW>
    err = abs(ex(i));     
end
if i==maxnoisesim
    [err,i] = min(abs(ex));
    data.noise   = noise{i};
    [~,~,~,dv] = model(params,data);
    x=[dv -2*S+3];
end

% bb = glmfit([dv -2*S+3],datx.actions==1,'binomial','constant','off');
% e1 = bb-[beta;c];

data.X       = struct('dQ',dv,'noisesim',i,'tolerr',tolerr,'maxnoisesim',maxnoisesim,'x',x,'err',err);

end

function [T,startpoint,endpoint,sdev,nS]=task_profile
T = 100;
nS = 2;

startpoint = -1;
endpoint   = +1;
sdev = 1;
end

function [R,walk]=randomwalk(T,startpoint,endpoint,s,nS)
if nargin<5, nS=1; end

N    = T/nS;
walk = nan(T,2);
r    = nan(T,2);
S    = nan(T,1);

if nS==2    
    S  = ones(T,1);
    tt = randperm(N);
    S(tt)  = 2;
else
    error('!');
end

for i=1:nS
    tt = (i-1)*N + (1:N);
    S(tt) = i;    
    steps = randn(N,2)*s;
    % ensure the sum is exactly zero
    steps = bsxfun(@minus,steps , mean(steps));
    % add in a bias to each step.
    ds = (endpoint - startpoint);
    steps = bsxfun(@plus,steps , [ds -ds]/N);
    walk(tt,:) = cumsum(steps);
    r(tt,:) = steps;
end
R = struct('S',S,'r',r);
end

function hierlap_fixedgroup(data,models,flap,config2,free_group)
kk = 1:length(models);

for k= kk
    [fdir,fname] = fileparts(flap{k});
    model = models{k};
    fhierlap = fullfile(fdir,sprintf('fghier%s.mat',fname));
    config2.free_group = free_group{k};
    config2.free_groupvar = free_group{k};
    if config2.save_prog, config2.fname_prog   = fullfile(fdir,sprintf('hier%s_prog_%0.4f.mat',fname,now)); end
    if ~exist(fhierlap,'file') && exist(flap{k},'file')
        if exist(flap{k},'file')
            cbm = cbm_hierlap(data, model, flap{k}, [], config2); %#ok<NASGU>
            save(fhierlap,'cbm');    
        else
            pause(90);
        end
    end
end

end

function hierlap_fixedmean(data,models,flap,config2,free_group)
kk = 1:length(models);

for k= kk
    [fdir,fname] = fileparts(flap{k});
    model = models{k};
    fhierlap = fullfile(fdir,sprintf('fmhier%s.mat',fname));
    config2.free_group = free_group{k};
    config2.free_groupvar = ones(size(free_group{k}));
    if config2.save_prog, config2.fname_prog   = fullfile(fdir,sprintf('hier%s_prog_%0.4f.mat',fname,now)); end
    
%     if exist(fhierlap,'file')
%         delete(fhierlap);
%     end
    
    if ~exist(fhierlap,'file') && exist(flap{k},'file')
        if exist(flap{k},'file')
            cbm = cbm_hierlap(data, model, flap{k}, [], config2); %#ok<NASGU>
            save(fhierlap,'cbm');
        else
            pause(90);
        end
    end
end

end