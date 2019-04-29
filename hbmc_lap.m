function hbmc_lap(data,models,flap,init,v0,config1,config2,nn,m0)
if nargin<9, m0 = repmat({0},1,length(models)); end
    
N  = length(data);
kk = 1:length(models);
hierfit = 1;
if nargin>7, hierfit = 0; end

for k= kk
    [fdir,fname] = fileparts(flap{k});
    model = models{k};
    config1.save_data  = 0; config1.save_prog  = 0;  
    
    if ~exist(flap{k},'file')
        for n=nn
            d     = length(init{k});            
            mm0   = m0{k}.*ones(d,1);
            prior = struct('mean',mm0,'variance',v0);        

            fdirlaps     = fullfile(fdir,'laps'); makedir(fdirlaps);
            flapn        = fullfile(fdirlaps,sprintf('%s_%03d.mat',fname,n));
            if ~exist(flapn,'file')
                cbm   = cbm_lap(data(n), model, prior, [], config1);     %#ok<NASGU>
                save(flapn,'cbm');
            end
        end    
        [fnames,okN] = getfileordered(fdirlaps,sprintf('%s_%s.mat',fname,'%03d'),1:N);
        if okN
            cbm = cbm_lap_aggregate(fnames); %#ok<NASGU>
            save(flap{k},'cbm');
        end
    end
    if hierfit
        fhierlap = fullfile(fdir,sprintf('hier%s.mat',fname));
        if config2.save_prog, config2.fname_prog   = fullfile(fdir,sprintf('hier%s_prog_%0.4f.mat',fname,now)); end
        if ~exist(fhierlap,'file') && exist(flap{k},'file')
            if exist(flap{k},'file')
                cbm = cbm_hierlap(data, model, flap{k}, [], config2); %#ok<NASGU>
                save(fhierlap,'cbm');    
            else
                pause(90);
            end
        end
    end    
end


end
