function [data,models,flap,init,v0,pconfig,fsave,flog,fname,config] = sim_sim_dataload(simdir)
    data   = fullfile(simdir,'data.mat');                

    fdata  = load(data);
    data   = fdata.data;
    models = fdata.models;
    init   = fdata.init;

    K      = length(models);
    flap   = cell(K,1);
    mnames = cell(1,K);
    pconfig(1:K,1) = deal(cbm_config(length(init{1}),'hierlap',[]));           
    config = struct('verbose',0,'save_data',0);
    for k=1:K
        flap{k} = fullfile(simdir,sprintf('lap_model%d.mat',k));
        pconfig(k) = cbm_config(length(init{k}),'hierlap',config);

        if ishandle(models{k})
            mnames{k}  = func2str(models{k});
        elseif  ischar(models{k})
            mnames{k}  = (models{k});
        end
    end

    fsave   = [];
    flog    = fullfile(simdir,'hbmc.log');
    fname   = fullfile(simdir,'hbmc.mat');

    v0 = prior_v0(0.05);
end

function v0 = prior_v0(precision)
x         = log(1/precision-1);
v0        = fzero(@(v)log(normpdf(x,0,v)/normpdf(0,0,v))-log(.5),3)^2;   
end