function [mx,ex,mp,ep,methods] = sim_stat_C3(fsimfits)
nf = size(fsimfits,1);
mx = cell(1,nf);
ex = cell(1,nf);
mp = cell(1,nf);
ep = cell(1,nf);

for i=1:nf
    [mx{i},ex{i},mp{i},ep{i},methods] = stat(fsimfits(i,:));
end

end

function [mx,ex,mp,ep,methods] = stat(fsimfits)
nv = length(fsimfits);
exlap = nan(2,nv);
mxlap = nan(1,nv);
eplap = nan(2,nv);
mplap = nan(1,nv);

for j=1:nv
    simfit  = load(fsimfits{j});
    config  = simfit.config;
    fit     = simfit.fit;

    [mxl,plap] = stat_run(fit,config);
    
    exlap(:,j) = mxl([1 3]);
    mxlap(j)  = mxl(2);
    
    eplap(:,j) = plap([1 3]);
    mplap(j)  = plap(2);   
end

mx = mxlap;
ex = exlap;

mp = mplap;
ep = eplap;

methods = {'NHI'};
end

% function [mhbi,mlap,phbi,plap,mshbi,mslap] = stat_C4(fit,config)
function [mlap,plap] = stat_run(fit,config)
sigmoid = @(x)1./(1+exp(-x));
kref = 1;
iref = 1;
nsim = 100;

for i=1:nsim
    z = fit(i).sim.z(kref,:)==1;
    h = fit(i).sim.h{kref}(iref,z);
    h = sigmoid(h);
    
    a = fit(i).lap(kref).parameters(z,iref)';
    a = sigmoid(a);        
    elap(i,:) = mean(abs(a-h));     %#ok<AGROW>
    
    pxl(i,:) = fit(i).rfxlap.pxp(kref); %#ok<AGROW>
    
    [~,ib] = max(config.Nbar);    
end
mlap = quantile(elap',[.25 .5 .75]);

pxlb = (ib==kref).*pxl + (ib~=kref).*(1-pxl);

plap = quantile(pxlb',[.25 .5 .75]);
end