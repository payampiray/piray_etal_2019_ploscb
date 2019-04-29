function [pxpm,pxps,Nbarm,Nbars,ngm,egm,rmm,rms,mdx,sdx,corrBMS,methods] = sim_stat_A1(fit,config)

N       = config.N;
K       = length(config.Nbar);
normx   = config.normx;
nsim    = length(fit);
for n=1:nsim
    [dfxhbi(n,:),dfxhier(n,:),dfxlap(n,:)] = sim_statfx_error(fit(n),normx);
    pxp(n,:)    = fit(n).hbi.protected_exceedance_prob;
    Nbar(n,:)   = fit(n).hbi.model_frequency/N;

    sz     = fit(n).sim.z;
    [~,iz] = max(sz,[],1);        

    r      = fit(n).hbi.responsibility';
    [~,ix] = max(r,[],1);
    ixgood = find(ix==iz);
    ixbad  = find(ix~=iz);
    indgood= sub2ind([K,N],ix(ix==iz),ixgood);
    indbad = sub2ind([K,N],ix(ix~=iz),ixbad);

    ngood(n,:)   = mean(ix==iz);
    rm(n,:)      = [median(r(indgood),2) median(r(indbad),2)];

    %----
    im = find(config.Nbar>(N/2));    
    if isempty(im), im = 2; end
    
    pxplap = fit(n).rfxlap.pxp;
    [~,ilap] = max(pxplap);
    bms(n,1) = sum(im==ilap');
    
    Fhier = [];
    for k=1:K
        Fhier = [Fhier fit(n).hier.log_evidence];
    end
    [~,ihier] = max(Fhier);
    bms(n,2) = sum(im==ihier);
    
    [~,ihbi] = max(fit(n).hbi.protected_exceedance_prob);
    bms(n,3) = sum(im==ihbi);
end

%------------
% stats
pxpm  = mean(pxp);
pxps  = serr(pxp);

Nbarm = mean(Nbar);
Nbars = serr(Nbar);
rmm   = nanmean(rm);
rms   = nanserr(rm);
ngm   = mean(ngood);
egm   = serr(ngood);

corrBMS = sum(bms,1)/nsim*100;

for k=1:size(dfxhbi,2)
    dx_lap = cell2mat(dfxlap(:,k));
    dx_hier = cell2mat(dfxhier(:,k));
    dx_hbi = cell2mat(dfxhbi(:,k));
    
    mdx_lap = mean(dx_lap);
    mdx_hier = mean(dx_hier);
    mdx_hbi = mean(dx_hbi);
    
    sdx_lap = serr(dx_lap);
    sdx_hier = serr(dx_hier);        
    sdx_hbi = serr(dx_hbi);    
    
    mdx{k} = [mdx_lap; mdx_hier; mdx_hbi];
    sdx{k} = [sdx_lap; sdx_hier; sdx_hbi];
end
methods = {'NHI','HPE','HBI'};
end
