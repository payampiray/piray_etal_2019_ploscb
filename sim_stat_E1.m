function [mpxp, elpxp, ehpxp, mnbar, elnbar, ehnbar, ms95, ms50, mdx, eldx, ehdx, methods_rfx, X, Y, auc, methods] = sim_stat_E1(fsimfits)
kref = 1;
nokref = 3-kref;

N       = size(fsimfits,2);
mpxp    = nan(2,N);
elpxp   = nan(2,N);
ehpxp   = nan(2,N);
mnbar   = nan(2,N);
elnbar  = nan(2,N);
ehnbar  = nan(2,N);
ms95    = nan(2,N);
ms50    = nan(2,N);

mdx     = nan(3,N);
eldx    = nan(3,N);
ehdx    = nan(3,N);
X       = cell(3,N);
Y       = cell(3,N);
auc     = nan(3,N);


for i=1:N
    zsim  = cell(1,2);
    phbi  = cell(1,2);
    plap  = cell(1,2);
    phier = cell(1,2);
    nhbi  = cell(1,2);
    nlap  = cell(1,2);    
    phbib = cell(1,2);
    plapb = cell(1,2);    
    
    dxhbi = cell(1,2);
    dxlap = cell(1,2);
    dxhier= cell(1,2);
    
    
    
    for j=1:size(fsimfits,1)
        simfit  = load(fsimfits{j,i});
        config  = simfit.config;
        fit     = simfit.fit;

        [~,kb] = max(config.Nbar); 
        zsim{j} = (kref==kb)+zeros(1,length(fit));                
        
        for n=1:length(fit)
            
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
            [~,~,~,dxhbi{j}(n),dxhier{j}(n),dxlap{j}(n)] = sim_statfx_error(fit(n));
        end
        
        % PXP of the best model
        phbib{j} = (kb==kref)*phbi{j} + (kb~=kref)*(1-phbi{j});
        plapb{j} = (kb==kref)*plap{j} + (kb~=kref)*(1-plap{j});
        
    end
    z   = cell2mat(zsim);    
    

    % note: methods_rfx = {'NHI','HBI'};
    nb1  = cell2mat(nlap);
    nb3  = cell2mat(nhbi);
    nb   = [nb1; nb3];
    
    pb1 = cell2mat(plapb);
    pb3 = cell2mat(phbib);    
    pb  = [pb1; pb3];        
          
    for m=1:2
        mpxp(m,i)  = median(pb(m,:) );
        elpxp(m,i) = quantile(pb(m,:),.25);
        ehpxp(m,i) = quantile(pb(m,:),.75);

        mnbar(m,i)  = median(nb(m,:) );
        elnbar(m,i) = quantile(nb(m,:),.25);
        ehnbar(m,i) = quantile(nb(m,:),.75);        
    end
    ms95(:,i) = mean(pb>.95,2)*100;
    ms50(:,i) = mean(pb>.5,2)*100;
        
    % note: methods_rfx = {'NHI','HPE','HBI'};
    p1  = cell2mat(plap);
    p2  = cell2mat(phier);
    p3  = cell2mat(phbi);
    p   = [p1; p2; p3];    
    
    d1  = cell2mat(dxlap);
    d2  = cell2mat(dxhier);
    d3  = cell2mat(dxhbi);
    dx  = [d1; d2; d3];    
    
    for m=1:3
        [X{m,i},Y{m,i},~,auc(m,i)] = perfcurve(z,p(m,:),1);
% %         [prec, tpr, fpr, thresh] = prec_rec(p(m,:), z);
        
        mdx(m,i)  = median(dx(m,:) );
        eldx(m,i) = quantile(dx(m,:),.25);
        ehdx(m,i) = quantile(dx(m,:),.75);        
    end
end

methods_rfx = {'NHI','HBI'};
methods = {'NHI','HPE','HBI'};
end   
