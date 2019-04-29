function emp_PD(gname)

ii =input(1);

pipedir    = getdefaults('pipedir');
thisdir    = fullfile(pipedir,'emp_PD');
dataset    = fullfile(thisdir,'data.mat');
dataset    = load(dataset);

switch gname
    case 'healthy'
        iig = 1:20;
    case 'ON'
        iig = 46:76;
    case 'HON'
        iig = [1:20 46:76];
    otherwise
        error('!');        
end

nn = 1:length(iig);
mm = [1 2 3];

mnames     = {'WSLS','RL','Dual-\alpha RL'};
models     = {@wsls,@qlearning,@q2learning};
d          = [1 3 4];


%------------------
data       = dataset.subs(iig)';

mnames     = mnames(mm);
models     = models(mm);
d          = d(mm);

fitdir     = fullfile(thisdir,gname); makedir(fitdir);

normx  = cell(1,length(d));
pnames = cell(1,length(d));
modelstrs = cell(1,length(d));
for i=1:length(d)
    [~,normx{i},pnames{i}] = models{i}(rand(1,d(i)), data{1});
    modelstrs{i} = func2str(models{i});
end

config  = struct('models',{modelstrs},'mnames',{mnames},'pnames',{pnames},'normx',{normx}); %#ok<NASGU>
fconfig = fullfile(fitdir,'config.mat');
if ~exist(fconfig,'file')
    save(fconfig,'config' );
end
%--------------------

v0         = 6.25;
alg        = 'hierlap';
config     = struct('verbose',1,'functionname','','save_prog',1);
K = length(d);
init = cell(K,1);
flap = cell(K,1);        
pconfig(1:K,1) = deal(cbm_config(length(init{1}),alg,[]));           
for k=1:K
    init{k} = zeros(1,d(k));
    flap{k} = fullfile(fitdir,sprintf('lap_model%d.mat',k));
    pconfig(k) = cbm_config(length(init{k}),alg,config);
end

hbmcdir = fitdir;
fsave   = [];
flog    = fullfile(hbmcdir,'hbmc.log');
fname   = fullfile(hbmcdir,'hbmc.mat');
fname0  = fullfile(hbmcdir,'hbmc0.mat');
flog0   = fullfile(hbmcdir,'hbmc0.log');

switch ii
    case 0
        hbmc_lap(data,models,flap,init,v0,config,[],nn);
    case 1
        if ~exist(fname,'file')
            hbmc(data,models,flap,fname,pconfig,fsave,flog);
        end
    case 2
        if ~exist(fname0,'file')
            hbmc(data,models,flap,fname0,pconfig,fsave,flog0,1);
        end
        fhbi = fullfile(hbmcdir,'hbi.mat');
        if ~exist(fhbi,'file')
            cbm  = hbmc_exceedance(fname,fname0); %#ok<NASGU>
            save(fhbi,'cbm');
        end
    otherwise
        error('!');
end
end

function [loglik,normx,pnames] = q2learning(params, data)
fxu    = @(t)(1./(1+exp(-t)));
beta   = exp(params(1));
alpha1 = fxu(params(2));
alpha2 = fxu(params(3));
pers   = params(4);

normx  = {@exp,@(x)(1./(1+exp(-x))),@(x)(1./(1+exp(-x))),@(x)x};
pnames = {'\beta','\alpha^+','\alpha^-','\it p'};

% extracting information
data        = data.mat;
state       = data(:,2);
reward      = data(:,4);
punishment  = data(:,5);
action      = data(:,6);
correctness = data(:,7);

outcome = correctness.*reward - (1-correctness).*punishment;

qvalue  = zeros(max(state),max(action));

qpers = ones(max(state),max(action));
p     = nan(size(action,1),1);

for t=1:size(action,1)
    
    % computing nll
    a = action(t);
    s = state(t);
    o = outcome(t);
    
    dq = beta*(qvalue(s,1)-qvalue(s,2))+pers*(qpers(s,1)-qpers(s,2));
    q  = 1./(1+exp(-dq));
    p(t) = q*(a==1)+ (1-q)*(a==2);
        
    % decaying perseveration
    qpers(s,:) = 0; qpers(s,a)=1;
    
    % learning
    delta      = o - qvalue(s,a);
    alpha      = alpha1*(delta>=0)+alpha2*(delta<0);
    qvalue(s,a)= qvalue(s,a)  + alpha*delta;
end

loglik = sum(log(p+eps));

end

function [loglik,normx,pnames] = qlearning(params, data)
px     = params;
px(3)  = px(2);
px(4)  = params(3);
[loglik,normx,pnames] = q2learning(px, data);
normx  = normx([1 2 4]);
pnames = pnames([1 2 4]);
end

function [loglik,normx,pnames] = wsls(params, data)
px     = params;
px(2)  = inf;
px(3)  = inf;
px(4)  = 0;
[loglik,normx] = q2learning(px, data);
normx = normx(1);
pnames = {'\beta'};
end