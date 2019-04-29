function [r,Nbar,thetabar,Sdiag,bound] = hbmc_qHZ(qmutau,qm,cm)
qmlimInf  = qm.limInf;

logf      = cm.logf;
theta     = cm.theta;
logdetA   = cm.logdetA;
Ainvdiag  = cm.Ainvdiag;

[K,N]    = size(logf);
r        = nan(K,N);
thetabar = cell(K,1);
Sdiag    = cell(K,1);
Nbar     = nan(K,1);
ElogpH   = nan(K,1);
ElogpZ   = nan(K,1);
ElogpX   = nan(K,1);
ElogqH   = nan(K,1);
ElogqZ   = nan(K,1);

D        = arrayfun(@(k)length(qmutau(k).a),(1:K)');
ElogdetT = arrayfun(@(k)sum(qmutau(k).Elogtau),(1:K)');
logdetET = arrayfun(@(k)sum(log(qmutau(k).Etau)),(1:K)');
beta     = arrayfun(@(k)qmutau(k).beta,(1:K)');

lambda   = .5*ElogdetT -.5*logdetET -.5*D./beta;

logrho   = logf -.5*logdetA;
logrho   = bsxfun(@plus,logrho,.5*D*log(2*pi) +lambda +qm.Elogm);

if qmlimInf, r = ones(K,N)/K; end

logeps      = exp(log1p(-1+eps));
for k=1:K
    if ~qmlimInf
    rarg        = bsxfun(@minus,logrho,logrho(k,:));
    r(k,:)      = 1./sum(exp(rarg),1);    
    end
    
    Nk          = sum(r(k,:));
    Dk          = D(k);
    
    Nbar(k)     = Nk;    
    thetabar{k} = sum(bsxfun(@times,theta{k},r(k,:)),2)/Nk;
    Sdiag{k}    = sum(bsxfun(@times,theta{k}.^2+Ainvdiag{k} ,r(k,:)),2)/Nk -thetabar{k}.^2;
    
    ElogpH(k)   = +.5*Nk*ElogdetT(k) -.5*Nk*Dk*log(2*pi) -.5*Nk*Dk/qmutau(k).beta +...
                  -.5*Nk*sum(qmutau(k).Etau.*( Sdiag{k}+(thetabar{k}-qmutau(k).a).^2 ) );
    ElogpZ(k)   = Nk*qm.Elogm(k);
    
    ElogpXH     = sum( r(k,:).*(logf(k,:) -.5*Dk + lambda(k)) );
    ElogpX(k)   = ElogpXH - ElogpH(k);
    
    rlogr       = r(k,:).*(log1p(-1+r(k,:))); % =r(k,:).*log(r(k,:))
    rlogr(r(k,:)<logeps) = 0;
    ElogqH(k)   = sum( r(k,:).*(-Dk/2-Dk/2*log(2*pi)+.5*logdetA(k,:) ) );
    ElogqZ(k)   = sum( rlogr);  
end


bound = struct('module','qHZ','ElogpX',ElogpX,'ElogpH',ElogpH,'ElogpZ',ElogpZ,'ElogqH',ElogqH,'ElogqZ',ElogqZ);
end