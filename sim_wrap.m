function sim_wrap(nsim,simcat,simstr,modelnames,mnames,pnames,normx)

tempdir  = getdefaults('tempdir');

pipedir  = getdefaults('pipedir');
makedir(fullfile(pipedir,simcat));
fsimfit  = fullfile(pipedir,simcat,sprintf('fit_%s.mat',simstr));

fsimfullfit  = fullfile(pipedir,'..','fullfits',simcat); makedir(fsimfullfit);
fsimfullfit  = fullfile(fsimfullfit,sprintf('full_%s.mat',simstr));

exfit = exist(fsimfit,'file');

if ~exfit
    simcatdir = fullfile(tempdir,simcat,simstr);        
    ok = zeros(1,nsim);
    while ~all(ok)            
        try %#ok<TRYNC>
            for j=1:nsim
                simdir    = fullfile(simcatdir,sprintf('sim%02d',j));
                ok(j)=exist(fullfile(simdir,'hbmc0.mat'),'file');
            end
%         catch
%             error('not all files exist');
        end
    end
    for j=1:nsim
        simdir    = fullfile(simcatdir,sprintf('sim%02d',j ));
        [config,fit(j),fullfit(j)] = wrapit(simdir); %#ok<AGROW,NASGU>
    end
    config.modelnames = modelnames;
    config.mnames     = mnames;
    config.pnames     = pnames;
    config.normx      = normx;    
    
    save(fsimfit,'config','fit');
    save(fsimfullfit,'config','fit','fullfit');
end
end

function [config,fit,fullfit] = wrapit(simdir)

fhbi    = fullfile(simdir,'hbi.mat');
if ~exist(fhbi,'file')
fname   = fullfile(simdir,'hbmc.mat');
math    = load(fname); math = math.math;
cbm     = math(end); %#ok<NASGU>
save(fhbi,'cbm');
end
fhbi    = load(fhbi); cbm = fhbi.cbm;

try
fname   = fullfile(simdir,'hbmc0.mat');
math0   = load(fname); math0 = math0.math;
cbm0    = math0(end);
[cbm]   = hbmc_exceedance(cbm,cbm0);
catch
error('No hbmc0');
end

cbm = hbmc_errorbar(cbm);

cbmoutput = hbmc_output(cbm);
hbi = cbmoutput;
hbi.logf = cbm.cm.logf';
hbi.math = cbm;

Nbar    = cbmoutput.model_frequency;
K       = size(Nbar,2);

fg = 0;
fm = 0;

% hierlap
for k=1:K

    flap = fullfile(simdir,sprintf('lap_model%d.mat',k));
    flap  = load(flap); flap = flap.cbm;    
    
    lmelap(:,k) = flap.output.log_evidence; %#ok<AGROW>
    
    lapfull(k) = flap; %#ok<AGROW>
    lap(k) = flap.output; %#ok<AGROW>
        
    fhier = fullfile(simdir,sprintf('hierlap_model%d.mat',k));
    fhier  = load(fhier); fhier = fhier.cbm;
    hierfull(k) = fhier; %#ok<AGROW>
    
    hierk = fhier.output;
    hierk.logf = fhier.math.quad_apx.logf';
    hier(k) = hierk; %#ok<AGROW>
    
    fhier = fullfile(simdir,sprintf('fghierlap_model%d.mat',k));
    if exist(fhier,'file')
        fg = 1;
        fhier  = load(fhier); fhier = fhier.cbm;
        hierk = fhier.output;
        hierk.logf = fhier.math.quad_apx.logf';
        fghier(k) = hierk; %#ok<AGROW>
        
        fghierfull(k) = fhier; %#ok<AGROW>
    end
    
    fhier = fullfile(simdir,sprintf('fmhierlap_model%d.mat',k));
    if exist(fhier,'file')
        fm = 1;
        fhier  = load(fhier); fhier = fhier.cbm;
        hierk = fhier.output;
        hierk.logf = fhier.math.quad_apx.logf';
        fmhier(k) = hierk; %#ok<AGROW>
        
        fmhierfull(k) = fhier; %#ok<AGROW>
    end    
end

[alpha,exp_r,xp,pxp,bor,g] = cbm_spm_BMS(lmelap, 1e6, 0);
rfxlap = struct('alpha',alpha,'exp_r',exp_r,'xp',xp,'pxp',pxp,'bor',bor,'g',g,'lme',lmelap);

fsim = fullfile(simdir,'sim.mat');
fsim = load(fsim); 
fsim = fsim.sim;
sim  = struct('h',{fsim.h},'z',fsim.z);

config = fsim.config;

fdata = fullfile(simdir,'data.mat');
fdata = load(fdata); 
data  = fdata.data;


fit = struct('sim',sim,'lap',lap,'hier',hier,'hbi',hbi,'rfxlap',rfxlap);
if fg
    fit.fghier = fghier;
end
if fm
    fit.fmhier = fmhier;
end

if nargout>2
fullfit = struct('sim',sim,'data',{data},'lap',lapfull,'hier',hierfull,'hbi',hbi);
    if fg
        fullfit.fghier = fghierfull;
    end
    if fm
        fullfit.fmhier = fmhierfull;
    end    
end

end