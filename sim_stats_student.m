function [ahbi,sehbi,Nbar,xmsim,sesim,dfsim,xmlap,selap,dflap,Fhier,xmhier,sehier,dfhier]=sim_stats_student(fsimfit,kref,ipref)

sumdir = getdefaults('sumdir');


[simcat,fname] = fileparts(fsimfit);
[~,simcat] = fileparts(simcat);
makedir(fullfile(sumdir,simcat));
fsum = fullfile(sumdir,simcat,sprintf('ttest_%s.mat',fname));

if ~exist(fsum,'file')
    simfit  = load(fsimfit);
    fit     = simfit.fit;
    zk      = fit(1).sim.z(kref,:)==1;
    nsim    = length(fit);

    xmsim  = nan(1,nsim);
    sesim  = nan(1,nsim);
    xmlap  = nan(1,nsim);
    selap  = nan(1,nsim);
    xmhier = nan(1,nsim);
    sehier = nan(1,nsim);    
    
    Fhier = nan(1,nsim);
    Nbar  = nan(1,nsim);
    ahbi  = nan(1,nsim);
    sehbi = nan(1,nsim);
    bor   = nan(1,nsim);


    for n=1:nsim                
        [xmsim(n),sesim(n),dfsim] = stats(fit(n).sim.h{kref}(ipref,zk)');
        [xmlap(n),selap(n),dflap] = stats(fit(n).lap(kref).parameters(:,ipref));

        bor(n) = fit(n).hbi.math.exceedance.bor;
        Nbar(n) = fit(n).hbi.model_frequency(kref);
        ahbi(n) = fit(n).hbi.group_mean{kref}(ipref);
        sehbi(n) = fit(n).hbi.group_errorbar{kref}(ipref);

        Fmain = fit(n).hier(kref).log_evidence;
        Fnull = fit(n).fmhier(kref).log_evidence;

        Fhier(n) = Fmain - Fnull;
        
        [xmhier(n),sehier(n),dfhier] = stats(fit(n).fmhier(kref).parameters(:,ipref));
    end

    st.sim = struct('xmsim',xmsim,'dfsim',dfsim,'sesim',sesim);
    st.lap = struct('xmlap',xmlap,'dflap',dflap,'selap',selap);
    st.hier = struct('Fhier',Fhier,'xmhier',xmhier,'dfhier',dfhier,'sehier',sehier);
    st.hbi = struct('ahbi',ahbi,'Nbar',Nbar,'sehbi',sehbi,'bor',bor); %#ok<STRNU>
    
    save(fsum,'-struct','st');
else
    sums  = load(fsum);

    xmsim = sums.sim.xmsim;
    sesim = sums.sim.sesim;
    dfsim = sums.sim.dfsim;    
    
    xmlap = sums.lap.xmlap;
    selap = sums.lap.selap;    
    dflap = sums.lap.dflap;   
    
    
    xmhier= sums.hier.xmhier;
    sehier= sums.hier.sehier;    
    dfhier= sums.hier.dfhier;
    Fhier = sums.hier.Fhier;
    
    ahbi  = sums.hbi.ahbi;
    sehbi = sums.hbi.sehbi;
    Nbar  = sums.hbi.Nbar;    
end

end

function [xmean,ser,df]=stats(x)
dim = 1;
nans = isnan(x);
if any(nans(:))
    samplesize = sum(~nans,dim);
else
    samplesize = size(x,dim); % a scalar, => a scalar call to tinv
end
df = max(samplesize - ones('like',x), 0); % make sure df is the same type as X
xmean = nanmean(x,dim);
sdpop = nanstd(x,[],dim);
sqrtn = sqrt(samplesize);

ser = sdpop ./ sqrtn;
end