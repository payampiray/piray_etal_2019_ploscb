function [qm,bound] = hbmc_qm(pm,Nbar)

limInf     = pm.limInf;
logC0      = pm.logC;
alpha0     = pm.alpha;

alpha      = alpha0 + Nbar;
alpha_star = sum(alpha);
Elogm      = psi(alpha) - psi(alpha_star);
loggamma   = gammaln(alpha);
logC       = gammaln(alpha_star)-sum(loggamma);

Elogpm     = logC0 + sum((alpha0-1).*Elogm);
Elogqm     = logC  + sum((alpha-1).*Elogm);
ElogpZ     = Nbar.*Elogm;

if limInf
    K        = length(alpha);
    pm.alpha = inf(K,1);
    Elogm    = log(ones(K,1)/K);
    logC     = inf;
    Elogpm   = nan;
    Elogqm   = nan;
    ElogpZ   = Nbar.*Elogm;
end

qm         = struct('name','Dirichlet','limInf',limInf,'alpha',alpha,'Elogm',Elogm,'logC',logC);
bound      = struct('module','qm','ElogpZ',ElogpZ,'Elogpm',Elogpm,'Elogqm',Elogqm);

end