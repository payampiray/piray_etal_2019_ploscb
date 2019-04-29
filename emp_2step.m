function emp_2step(input)
ii =input(1);
nn = 0;
if length(input)==2, nn=input(2); end
run(ii,nn);
end

function run(ii,nn)
pipedir    = getdefaults('pipedir');
fitdir     = fullfile(pipedir,'emp_2step');
[data,models,flap,init,v0,pconfig,fsave,flog,fname,configlap,config] = dataload(fitdir); %#ok<ASGLU>
% config.numinit_med = 30;
% config.numinit_up  = 30;
% config.discard_bad = 1;

fconfig = fullfile(fitdir,'config.mat');
if ~exist(fconfig,'file')
    save(fconfig,'config' );
end

switch ii
    case 0
        hbmc_lap(data,models,flap,init,v0,configlap,[],nn);
    case {1,2,3}
        models = models(ii);
        flap   = flap(ii);
        init   = init(ii);
        pconfig = pconfig(ii);
        hbmc_lap(data,models,flap,init,v0,[],pconfig);
    case 4
        hbmc(data,models,flap,fname,pconfig,fsave,flog);
    case 5 
        [fdir,fname0]=fileparts(fname); fname0=fullfile(fdir,sprintf('%s0.mat',fname0));
        [fdir,flog0]=fileparts(flog); flog0=fullfile(fdir,sprintf('%s0.mat',flog0));
        if ~exist(fname0,'file')
        hbmc(data,models,flap,fname0,pconfig,[],flog0,1);
        end
        
        fhbi = fullfile(fdir,'hbi.mat');
        if ~exist(fhbi,'file')
            cbm  = hbmc_exceedance(fname,fname0); %#ok<NASGU>
            save(fhbi,'cbm');
        end        
    otherwise
        error('!');
end
end

function [data,models,flap,init,v0,pconfig,fsave,flog,fname,configlap, config] = dataload(fitdir)

data = fullfile(fitdir,'data.mat');
data = load(data); data = data.data;
d = [7 4 6];

alg    = 'hierlap';
%         config = struct('verbose',0,'loop_runtime',55,'hessian','on','gradient','on','largescale','on');
% configlap = struct('verbose',1,'hessian','on','gradient','on','largescale','on','functionname','','save_prog',1);
configlap = struct('verbose',0);
mnames = {'Hybrid','MB','MF'};
% models = {@hbmc_2step_hybrid;@hbmc_2step_mb;@hbmc_2step_mf};
models = {@model_hybrid;@model_mb;@model_mf};

K = length(d);
init = cell(K,1);
flap = cell(K,1);        
pconfig(1:K,1) = deal(cbm_config(length(init{1}),alg,[]));           
for k=1:K
    init{k} = zeros(1,d(k));
    flap{k} = fullfile(fitdir,sprintf('lap_%s.mat',mnames{k}));
%             config.functionname = mnames{k};            
    pconfig(k) = cbm_config(length(init{k}),alg,configlap);
end     

hbmcdir = fitdir;
fsave   = [];
flog    = fullfile(hbmcdir,'hbmc.log');
fname   = fullfile(hbmcdir,'hbmc.mat');

v0      = 6.25;


fxnames    = {'\alpha_1', '\alpha_2', '\lambda', '\beta_1', '\it w', '\it p', '\beta_2'};
pnames     = { fxnames , fxnames([2 4 6 7]), fxnames([1 2 3 4 6 7]) };

normx1 = {@(x)1./(1+exp(-x)), @(x)1./(1+exp(-x)), @(x)1./(1+exp(-x)), @exp, @(x)1./(1+exp(-x)), @(x)x, @exp};
normx2 = normx1([2 4 6 7]);
normx3 = normx1([1 2 3 4 6 7]);

normx  = {normx1, normx2, normx3};

config  = struct('models',{models},'mnames',{mnames},'pnames',{pnames},'normx',{normx});

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
