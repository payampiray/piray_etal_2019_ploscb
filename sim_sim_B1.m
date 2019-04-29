function sim_sim_B1(ii,simcat,simstr,models)
    
for i=ii
    simrun(simcat,simstr,i,models);
end

end

function simrun(simcat,simstr,id,models)

simdir  = fullfile(getdefaults('tempdir'),simcat,simstr,sprintf('sim%02d',id));
fprintf('sim directory is %s\n',simdir);

for k=1:length(models)
    models{k} = str2func(models{k});
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

function [loglik,data,p] = model_kf1(params,dat)
gamma   = 0;
omega   = exp(params(1));
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
v = [1 1];
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
    
    vhat     = v(a)+gamma;
    kalman   = vhat./(vhat+omega);
    delta    = o - q(a)';    
    q(a)     = q(a) + (kalman.*delta)';
    v(a)     = (1-kalman).*vhat;
    
    actions(t) = a;
    outcome(t) = o;    
end

if simmode
    data = struct('actions',actions,'outcome',outcome,'simmed',dat);
end

loglik = sum(log(f+eps));
end
%--------------------------
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

function [R,data]=task_kf1(r,params)
gamma   = 0;
omega   = exp(params(1));

T = length(r);

x = [0 0];
X = nan(T,2);
v = [1 1];
for t=1:T
    delta = (r(t,:)-x);
%     delta2 = (r(t)-y);
%     alpha  = (delta2>=0)*alpha1 + (delta2<0)*alpha2;
    vhat   = v+gamma;
    kalman = vhat./(vhat+omega);
    x = x + kalman.*delta;    
    v = (1-kalman).*vhat;
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
