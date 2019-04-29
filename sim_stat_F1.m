function [mdgx,sdgx,pnames,methods]=sim_stat_F1(fsimfits,kref)

[nf, no]   = size(fsimfits);

simfit     = load(fsimfits{1});
config     = simfit.config;
D          = length(config.mu{kref});

mdgx = repmat({nan(3,no)},nf,D);
sdgx = repmat({nan(3,no)},nf,D);

normx     = config.normx;
pnames    = config.pnames{kref};

for f=1:nf
    for o=1:no
        fsimfit = fsimfits{f,o};
        simfit     = load(fsimfit);
        fit        = simfit.fit;        

        nsim    = length(fit);                
        for n=1:nsim            
            glap(n,:)  = mean(fit(n).lap(kref).parameters); %#ok<AGROW>
            ghier(n,:) = mean(fit(n).hier(kref).parameters); %#ok<AGROW>
            ghbi(n,:)  = fit(n).hbi.group_mean{kref}; %#ok<AGROW>
            gsim(n,:)  = config.mu{kref}'; %#ok<AGROW>        
        end
        
        %------
        % group parameters
        for d=1:D
            glap(:,d)  = normx{kref}{d}(glap(:,d));
            ghier(:,d)  = normx{kref}{d}(ghier(:,d));
            ghbi(:,d)  = normx{kref}{d}(ghbi(:,d));
            gsim(:,d)  = normx{kref}{d}(gsim(:,d));                        
        end     
        
        mdglap  = mean(abs(glap-gsim));
        edglap  = serr(abs(glap-gsim));
        mdghier = mean(abs(ghier-gsim));
        edghier = serr(abs(ghier-gsim));
        mdghbi  = mean(abs(ghbi-gsim));
        edghbi  = serr(abs(ghbi-gsim));
                
        
        mdxgfo    = [mdglap; mdghier; mdghbi];
        sdxgfo    = [edglap; edghier; edghbi];

%         mdxfo = [mdx_lap; mdx_hier; mdx_hbi];
%         sdxfo = [sdx_lap; sdx_hier; sdx_hbi];
        
        for d=1:D
%         mdx{f,d}(:,o) = mdxfo(:,d);
%         sdx{f,d}(:,o) = sdxfo(:,d);
        
        mdgx{f,d}(:,o) = mdxgfo(:,d);
        sdgx{f,d}(:,o) = sdxgfo(:,d);
        end        
    end    
end
methods = {'NHI','HPE','HBI'};    

end