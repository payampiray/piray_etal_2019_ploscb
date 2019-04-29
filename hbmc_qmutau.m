function [qmutau,pmutau,bound,optim] = hbmc_qmutau(pmutau,Nbar,thetabar,Sdiag,config)

b   = pmutau(1).beta;
v   = pmutau(1).nu;

%----------------------
K       = size(thetabar,1);
a       = cell(K,1);
beta    = cell(K,1);
sigma   = cell(K,1);
nu      = cell(K,1);
Etau    = cell(K,1);
Elogtau = cell(K,1);
logG    = cell(K,1);

for k=1:K
    a0k           = pmutau(k).a;
    beta0k        = b;
    nu0k          = v;
    sigma0k       = pmutau(k).sigma;
    
    Nk            = Nbar(k);
    
    beta{k}       = beta0k  + Nk;
    a{k}          = (beta0k*a0k + Nk*thetabar{k})/beta{k};
    nu{k}         = nu0k + .5*Nk;
    sigma{k}      = sigma0k + .5*( Nk*Sdiag{k} + Nk*beta0k/(Nk+beta0k)*(thetabar{k}-a0k).^2 );
    
    Elogtau{k}    = psi(nu{k})-log(sigma{k});
    Etau{k}       = nu{k}./sigma{k};
    logG{k}       = sum(-gammaln(nu{k}) + nu{k}*log(sigma{k}));
end
qmutau   = struct('name','GaussianGamma','a',a,'beta',beta,'sigma',sigma,'nu',nu,'Etau',Etau,'Elogtau',Elogtau,'logG',logG);

[~,bound] = bounding(qmutau,pmutau,Nbar,thetabar,Sdiag,b,v);
optim = struct('config',config);
end

function [L,bound] = bounding(qmutau,pmutau,Nbar,thetabar,Sdiag,b,w)
K          = size(thetabar,1);

ElogpH     = nan(K,1);
Elogpmu    = nan(K,1);
Elogqmu    = nan(K,1);
Elogptau   = nan(K,1);
Elogqtau   = nan(K,1);
for k=1:K
    a0k     = pmutau(k).a;
    beta0k  = b;
    nu0k = w;
    sigma0k = pmutau(k).sigma;
    
    logG0         = sum(-gammaln(nu0k) + nu0k.*log(sigma0k));
    
    Nk            = Nbar(k);
    Dk            = length(a0k);
    
    beta          = qmutau(k).beta;
    a             = qmutau(k).a;
    nu            = qmutau(k).nu;
    
    Elogtau       = qmutau(k).Elogtau;
    ElogdetT      = sum(Elogtau);
    Etau          = qmutau(k).Etau;
    ET            = diag(Etau);
    logG          = qmutau(k).logG;
    
    Elogpmu(k)    = -Dk/2*log(2*pi) +.5*Dk*log(beta0k) + .5*ElogdetT -.5*(a-a0k)'*(beta0k*ET)*(a-a0k) -Dk/2*beta0k/beta;
    Elogptau(k)   = +(nu0k-1)*ElogdetT -sum(sigma0k.*diag(ET)) + logG0;
    
    Elogqmu(k)    = -Dk/2*log(2*pi) +.5*Dk*log(beta) + .5*ElogdetT + -Dk/2;
    Elogqtau(k)   = +(nu-1)*ElogdetT - Dk*nu + logG;
    
    ElogpH(k)     = +.5*Nk*ElogdetT -.5*Nk*Dk*log(2*pi) -.5*Nk*Dk/beta +...
                    -.5*sum(Etau.*( Nk*Sdiag{k}+Nk*(thetabar{k}-a).^2 ) );
end
L       = sum(ElogpH + Elogpmu + Elogptau - Elogqmu - Elogqtau);
bound   = struct('module','qmutau','ElogpH',ElogpH,'Elogpmu',Elogpmu,'Elogptau',Elogptau,'Elogqmu',Elogqmu,'Elogqtau',Elogqtau);
end
