function [ms, mdx, edx, methods] = sim_stat_D1(fsimfits)
kref = 1;
nokref = 3-kref;

nf       = size(fsimfits,2);
ms      = nan(3,nf);

mdx     = nan(3,nf);
edx     = nan(3,nf);
X       = cell(3,nf);
Y       = cell(3,nf);
auc     = nan(3,nf);


for i=1:nf
    zsim  = cell(1,2);
    phbi  = cell(1,2);
    plap  = cell(1,2);
    phier = cell(1,2);
    nhbi  = cell(1,2);
    nlap  = cell(1,2);    
    
    dxhbi = cell(1,2);
    dxlap = cell(1,2);
    dxhier= cell(1,2);
    
    
    
    for j=1:size(fsimfits,1)
        simfit  = load(fsimfits{j,i});
        config  = simfit.config;
        fit     = simfit.fit;

        N = config.N;
        K = config.K;
        [~,kb] = max(config.Nbar); 
        zsim{j} = (kref==kb)+zeros(1,length(fit));                
        
        nsim = length(fit);
        bms  = nan(nsim,1);
        
        for n=1:nsim
            
            % PXP of the ref model
            phbi{j}(n) = fit(n).hbi.protected_exceedance_prob(kref);
            plap{j}(n) = fit(n).rfxlap.pxp(kref);
            
            % Nbar of the best model
            nhbi{j}(n) = fit(n).hbi.model_frequency(kb)/sum(fit(n).hbi.model_frequency);
            nlap{j}(n) = fit(n).rfxlap.exp_r(kb);

            % evidence of the ref model
            Fk   = fit(n).hier(kref).log_evidence;
            Fnok = fit(n).hier(nokref).log_evidence;
            Fhier = Fk - Fnok;    
            phier{j}(n) = 1./(1+exp(-Fhier));
            
            % mean error
            [d1,d2,d3,dxhbi{j}(n),dxhier{j}(n),dxlap{j}(n)] = sim_statfx_error(fit(n));
            
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
                
    end
    
    ms(:,i) = sum(bms,1)'/nsim*100;
    
    d1  = cell2mat(dxlap);
    d2  = cell2mat(dxhier);
    d3  = cell2mat(dxhbi);
    dx  = [d1; d2; d3];    
    
    for m=1:3
        mdx(m,i)  = mean(dx(m,:),2 );
        edx(m,i)  = serr(dx(m,:),2 );
    end    
end

methods = {'NHI','HPE','HBI'};
end   
