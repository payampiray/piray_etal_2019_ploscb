function [corrBMS,methods] = sim_stat_B1(fsimfits)

corrBMS = nan(length(fsimfits),3);
for i=1:length(fsimfits)
    simfit = load(fsimfits{i});
    fit = simfit.fit;
    config = simfit.config;
    N       = config.N;
    K       = length(config.Nbar);
    nsim    = length(fit);
    for n=1:nsim
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
    corrBMS(i,:) = sum(bms,1)/nsim*100;
end
methods = {'NHI','HPE','HBI'};
end
