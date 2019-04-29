function output = hbmc_output(math)
cbm = math(end);

theta = cbm.cm.theta;
L = cbm.bound.bound.L;
Nbar = cbm.Nbar';
r = cbm.r';

K  = length(cbm.qmutau);
a  = cell(1,K);
se = cell(1,K);
nk = nan(1,K);
for k=1:K
    theta{k} = theta{k}';
    
    a{k}    = cbm.qmutau(k).a';
    se{k}   = cbm.qmutau(k).se'; 
    nk(k)   = cbm.qmutau(k).nk';
end

xp  = [];
pxp = [];
if isfield(cbm,'exceedance')
    xp = cbm.exceedance.xp;
    if isfield(cbm.exceedance,'pxp')
        pxp = cbm.exceedance.pxp;
    else
        pxp = nan(size(xp));
    end
end

output     = struct('parameters',{theta},'responsibility',r,...
                    'group_mean',{a},'group_errorbar',{se},'degrees_of_freedom',nk,...
                    'log_evidence',L,'model_frequency',Nbar,'exceedance_prob',xp,'protected_exceedance_prob',pxp);


end