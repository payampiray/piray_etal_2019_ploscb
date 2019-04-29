function hbmc_permute(data,models,flaps,fdir,pconfig,group,ii)

N = length(data);
ng = max(group);

K = length(models);

Ng = nan(1,ng);
for g=1:ng
    Ng(g) = sum(group==g);
end


for i=ii
    nn = randperm(N);
    npre = 0;    
    permgroup = nan(1,N);    
    for g=1:ng
        nng = nn( npre + (1:Ng(g)) ); npre = npre + Ng(g);
        
        permgroup(nng) = g;
    end    
    
    gflap = cell(K,1);
    for g=1:ng
        permutations = find(permgroup==g); %#ok<NASGU>
        gdata  = data(permgroup==g);
        
        
        for k=1:K
            cbm = cbm_lap_aggregate(flaps(permgroup==g,k)); %#ok<NASGU>
            gflap{k} = fullfile(fdir,sprintf('lap_model%d_perm%04d.mat',k,i));
            save(gflap{k},'cbm');        
        end
        
        fname = fullfile(fdir,sprintf('hbmc_perm%04d_g%d.mat',i,g));
        
        math = hbmc(gdata,models,gflap,[],pconfig,[],[]);
        cbm  = math(end); %#ok<NASGU>
        save(fname,'cbm','permutations');
    end
end



end