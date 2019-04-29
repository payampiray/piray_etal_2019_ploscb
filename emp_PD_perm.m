function emp_PD_perm(ii,gname,permname)

pipedir    = getdefaults('pipedir');
thisdir    = fullfile(pipedir,'emp_PD');
dataset    = fullfile(thisdir,'data.mat');
dataset    = load(dataset);

gnames     = {'healthy','ON'};
iig        = {1:20,46:76};


group = [];
nn = [];
ngindex = cell(1,length(gname));
for g=1:length(gname)
    ig         = iig{strcmp(gnames,gname{g})};
    nn         = [nn ig]; %#ok<AGROW>
    ng         = length(ig);
    group      = [group g*ones(1,ng)]; %#ok<AGROW>
    
    ngindex{g} = 1:length(ig);
end

% mnames     = {'WSLS','RL','Dual-\alpha RL'};
models     = {@wsls,@qlearning,@q2learning};
d          = [1 3 4];

data       = dataset.subs(nn)';

fitdir     = fullfile(thisdir,permname); makedir(fitdir);

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

N = length(nn);
flaps = cell(N,K);
ngindex = cell2mat(ngindex);

for n=1:N
    gn = group(n);
    for k=1:K
        flaps{n,k} = fullfile(thisdir,gname{gn},'laps',sprintf('lap_model%d_%03d.mat',k,ngindex(n) ) );
%         ex(n,k) = exist(flaps{n,k},'file');
    end
end

tempdir    = getdefaults('tempdir');
fdir = fullfile(tempdir,'emp_PD',permname);makedir(fdir);

hbmc_permute(data,models,flaps,fdir,pconfig,group,ii);

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
% qvalue(1:2,:) = +12.5;
% qvalue(3:4,:) = -12.5;

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