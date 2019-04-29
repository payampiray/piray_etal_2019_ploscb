function [mpxp,spxp,mnbar,snbar,mbms,mdmx,sdmx, mdx,sdx,methods] = sim_stat_G1(fsimfits)

for i=1:length(fsimfits)
    simfit = load(fsimfits{i});
    fit    = simfit.fit;
    config = simfit.config;
    
    N       = config.N;    
    K       = length(config.Nbar);
    [~,kref]= max(config.Nbar);
    normx   = config.normx;
    nsim    = length(fit);
    
    D = length(config.mu{kref});
    dxhbi = nan(nsim,D);
    dxhier = nan(nsim,D);
    dxlap = nan(nsim,D);
    for n=1:nsim
        [dfxhbi,dfxhier,dfxlap] = sim_statfx_error(fit(n),normx);
        dxhbi(n,:)  = dfxhbi{kref};
        dxhier(n,:) = dfxhier{kref};
        dxlap(n,:)  = dfxlap{kref};
        
        [~,~,~,dxmhbi(n,1),dxmhier(n,1),dxmlap(n,1)] = sim_statfx_error(fit(n));
        
        pxphbi(n,:)    = fit(n).hbi.protected_exceedance_prob;
        Nbarhbi(n,:)   = fit(n).hbi.model_frequency/N;
        %----
        im = find(config.Nbar>(N/K));    
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
    
    mdxhbi  = mean(dxhbi,1);
    mdxhier = mean(dxhier,1);
    mdxlap  = mean(dxlap,1);
    sdxhbi  = serr(dxhbi,1);
    sdxhier = serr(dxhier,1);
    sdxlap  = serr(dxlap,1);
    
    mdx{i}     = [mdxlap; mdxhier; mdxhbi];
    sdx{i}     = [sdxlap; sdxhier; sdxhbi];
    
    mdmx{i}    = [mean(dxmlap); mean(dxmhier); mean(dxmhbi)];
    sdmx{i}    = [serr(dxmlap); serr(dxmhier); serr(dxmhbi)];
    
    mbms{i}    = 100*mean(bms,1)';
    
    mpxp{i}    = mean(pxphbi,1);
    mnbar{i}   = mean(Nbarhbi,1);
    spxp{i}    = serr(pxphbi,1);
    snbar{i}   = serr(Nbarhbi,1);
    
    pnamesref(i)  = config.pnames(kref);
    mnamesref(i)  = config.mnames(kref);
    modelnamesref(i)  = config.modelnames(kref);
end
    
methods = {'NHI','HPE','HBI'};

mpxp = cell2mat(mpxp')';
spxp = cell2mat(spxp')';
mpxp(mpxp<.01)=0.01;

mnbar = cell2mat(mnbar')';
snbar = cell2mat(snbar')';

mbms = cell2mat(mbms);
mbms(mbms<1) = 1;

mdmx = cell2mat(mdmx);
sdmx = cell2mat(sdmx);
end