function [mpxp,spxp,mNbar,sNbar,auc,ms50,ms95,TR,ratings,xa,ya,auca,methods] = sim_stat_C2(fsimfits,kref)

zz = []; rr = []; gg = [];
for j=1:length(fsimfits)
    fsimfit = fsimfits{j};
    simfit  = load(fsimfit);
    fit     = simfit.fit;
    config  = simfit.config;
    
    [~,kb]  = max(config.Nbar);
    if diff(config.Nbar)==0, kb = 2; end
    
    z = [];
    r = [];
    g = [];

    for i=1:length(fit)
    %     zk = p;

        zi = fit(i).sim.z(kref,:)'==1;
        z = [z; zi];

        ri = fit(i).hbi.responsibility(:,kref);
        r = [r; ri];

        gi = fit(i).rfxlap.g(:,kref);
        g = [g; gi];

        
        pxp(i,1)    = fit(i).hbi.protected_exceedance_prob(kref);
        pxplap(i,1)  = fit(i).rfxlap.pxp(kref);
        
        pxpb(i,1)    = fit(i).hbi.protected_exceedance_prob(kb);
        pxplapb(i,1)  = fit(i).rfxlap.pxp(kb);        
        
        Nbar(i,1)   = fit(i).hbi.model_frequency(kref)/sum(fit(i).hbi.model_frequency);        
        Nbarlap(i,1) = fit(i).rfxlap.exp_r(kref);    
    end

    zz = [zz; z];
    rr = [rr; r];
    gg = [gg; g];
    
    [xr{j},yr{j},~,auchbi] = perfcurve(z,r,1);
    [xg{j},yg{j},~,auclap] = perfcurve(z,g,1);

    pxpm = mean(pxp);
    pxps = serr(pxp);

    Nbarm = mean(Nbar);
    Nbars = serr(Nbar);

    pxplapm = mean(pxplap);
    pxplaps = serr(pxplap);

    Nbarlapm = mean(Nbarlap);
    Nbarlaps = serr(Nbarlap);


    %------------------
    mpxp(:,j)  = [pxplapm; pxpm ];
    spxp(:,j)  = [pxplaps;pxps];
    mNbar(:,j) = [Nbarlapm; Nbarm];
    sNbar(:,j) = [Nbarlaps; Nbars];
    auc(:,j)   = [auclap; auchbi];

    ms50(:,j)  = [mean(pxplapb>.5);mean(pxpb>.5)]*100;
    ms75(:,j)  = [mean(pxplapb>.75); mean(pxpb>.75)]*100;
    ms95(:,j)  = [mean(pxplapb>.95); mean(pxpb>.95)]*100;
end

x = [xg; xr];
y = [yg; yr];

tr = 0.50;
TP = sum( (zz==1) & (gg>tr));
FN = sum( (zz==1) & (gg<tr));    
FP = sum( (zz==0) & (gg>tr));
TN = sum( (zz==0) & (gg<tr));
TRlap  = (TP+TN)/(TP+TN+FP+FN);
TPRlap = TP/(TP+FN);
TNRlap = TN/(TN+FP);
FPRlap = FP/(TN+FP);

TP = sum( (zz==1) & (rr>tr));
FN = sum( (zz==1) & (rr<tr));    
FP = sum( (zz==0) & (rr>tr));
TN = sum( (zz==0) & (rr<tr));
TRhbi  = (TP+TN)/(TP+TN+FP+FN);
TPRhbi = TP/(TP+FN);
TNRhbi = TN/(TN+FP);
FPRhbi = FP/(TN+FP);

% TR = [ [TPRlap TNRlap TRlap]; [TPRhbi TNRhbi TRhbi] ];
% ratings = {'True positive','True negative',sprintf('Accuracy')};

TR = [ [TPRlap FPRlap TRlap]; [TPRhbi FPRhbi TRhbi] ];
ratings = {'TP','FP',sprintf('Accuracy')};

[xalap,yalap,~,auclap] = perfcurve(zz,gg,1);
[xahbi,yahbi,~,auchbi] = perfcurve(zz,rr,1);

% % [prec, tpr, fpr, thresh] = prec_rec(gg, zz);
% % [prec, tpr, fpr, thresh] = prec_rec(rr, zz);

% % % sgg = sort(gg,'descend');
% % % sgg = [1;sgg];
% % % lfp = nan(size(xalap));
% % % ltp = nan(size(xalap));
% % % for i=1:length(gg)
% % %     [lfp(i),ltp(i)]=tfrate(zz,gg,sgg(i) );
% % % end

xa = {xalap; xahbi};
ya = {yalap; yahbi};
auca = [auclap; auchbi];

methods = {'NHI','HBI'};
end

function [FP,TP]=tfrate(zz,rr,tr)
TP = mean( (zz==1) & (rr>tr));
FP = mean( (zz==0) & (rr>tr));
end