function [cm,thetabar,Sdiag,pmutau,pm,bound,prog]= hbmc_init(flap,a0,b,v,sigma0,limInf)
if nargin<6, limInf = 0; end

bb        = struct('ElogpX',nan,...
                   'ElogpH',nan,'ElogpZ',nan,'Elogpmu',nan,'Elogptau',nan,'Elogpm',0,...
                   'ElogqH',nan,'ElogqZ',nan,'Elogqmu',nan,'Elogqtau',nan,'Elogqm',0,...
                   'pmlimInf',limInf,'lastmodule','','L',nan);


%------------ cm init               
% initializing of G and H parameters
K        = length(flap);
theta    = cell(K,1);
Ainvdiag = cell(K,1);
logdetA  = cell(K,1);
V0       = cell(K,1);
logf     = cell(K,1);
D        = nan(K,1);
for k    = 1:K
    cbm_map     = load(flap{k}); cbm_map     = cbm_map.cbm;
    [logf{k},theta{k},Ainvdiag{k},logdetA{k},V0{k}] = init_GH(cbm_map);    
    D(k) = size(theta{k},1);
end
logf    = cell2mat(logf);
logdetA = cell2mat(logdetA);
cm      = struct('logf',logf,'theta',{theta},'Ainvdiag',{Ainvdiag},'logdetA',logdetA);

N        = size(logf,2);
thetabar = cell(K,1);
Sdiag    = cell(K,1);
for k=1:K
    thetabar{k} = sum(theta{k},2)/N;    
    Sdiag{k}    = sum(theta{k}.^2+Ainvdiag{k},2)/N -thetabar{k}.^2;        
end

% p(mu,tau)
a       = cell(K,1);
beta    = cell(K,1);
sigma   = cell(K,1);
nu      = cell(K,1);
alpha0  = ones(K,1);
for k  = 1:K
    if iscell(a0)
        ak = a0{k};
    else
        ak = a0.*ones(D(k),1);
    end
    
    a{k}     = ak;
    beta{k}  = b;
    sigma{k} = sigma0*ones(size(a{k}));%V0{k};
    nu{k}    = v;
end
pmutau   = struct('name','GaussianGamma','a',a,'beta',beta,'nu',nu,'sigma',sigma);
pmutau   = cal_GaussianGamma(pmutau);
pm       = struct('name','Dirichlet','limInf',limInf,'alpha',alpha0);
pm       = cal_Dirichlet(pm);

bound     = struct('bound',bb,'qHZ',[],'qmutau',[],'qm',[]);
prog      = struct('L',bound.bound.L,'r',ones(K,N)/K,'alpha',pm.alpha,'x',nan);
end

function [logf,theta,Ainvdiag,logdetA,V0] = init_GH(cbm)
N        = length(cbm.math.lme);
logf     = cbm.math.loglik;
% % logrho   = cbm.math.lme;
theta    = cell2mat(cbm.math.theta);
Ainvdiag = cell2mat(arrayfun(@(n)diag(cbm.math.T{n}^-1),1:N,'UniformOutput',0));
logdetA  = cell2mat(arrayfun(@(n)2*sum(log(diag(chol(cbm.math.T{n})))),1:N,'UniformOutput',0));
V0       = diag(cbm.input.prior.precision).^-1;
end

function [qmutau] = cal_GaussianGamma(qmutau)
K  = length(qmutau);
for k=1:K
    a         = qmutau(k).a;
    beta      = qmutau(k).beta;    
    nu        = qmutau(k).nu;
    sigma     = qmutau(k).sigma;    
    
    
    % B.30    
    Elogtau   = psi(nu)-log(sigma);
    
    % B.27
    Etau      = nu./sigma; 
    
    % logG
    logG      = sum(-gammaln(nu) + nu.*log(sigma));
    
    % update
    qmutau(k).a         = a;
    qmutau(k).beta      = beta;    
    qmutau(k).Etau      = Etau;
    qmutau(k).Elogtau   = Elogtau;
    qmutau(k).logG      = logG;
end

end

function [pm]     = cal_Dirichlet(pm)
% Bishop, B.21 for Dir(alpha)

alpha      = pm.alpha;
alpha_star = sum(alpha);
Elogm      = psi(alpha) - psi(alpha_star);

loggamma1  = gammaln(alpha);
logC       = gammaln(alpha_star)-sum(loggamma1);

if pm.limInf
    pm.alpha = inf(size(alpha));
    Elogm = inf(size(alpha));
    logC  = 0;
end

pm.Elogm   = Elogm;
pm.logC    = logC;

end