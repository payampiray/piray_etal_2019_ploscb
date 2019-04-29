function [pxp,Nbar,mnames,mx,ex,pnames,L,alpha]= emp_stat(fnames,fconfigs)

nf      = length(fnames);
pxp     = cell(nf,1);
Nbar    = cell(nf,1);
mnames  = cell(nf,1);
pnames  = cell(nf,1);
mx      = cell(nf,1);
ex      = cell(nf,1);

L = nan(nf,1);
for i=1:length(fnames)
    fname = fnames{i};
    fconf = fconfigs{i};
    fcbm  = load(fname); cbm = fcbm.cbm;
    fconfig = load(fconf); config = fconfig.config;


    Nbari   = cbm.Nbar'/sum(cbm.Nbar);
    pxpi    = cbm.exceedance.pxp;

    [~,im]  = max(Nbari);
    cbm     = hbmc_errorbar(cbm);
    qmutau  = cbm.qmutau(im);
    ahbi    = qmutau.a';
    sehbi   = qmutau.se';
    
    L(i) = cbm.bound.bound.L;
    alpha{i} = cbm.qm.alpha;

    mnamesi = config.mnames;
    pnamesi = config.pnames{im};
    normx   = config.normx{im};
    D       = length(ahbi);
    mxi     = nan(1,D);
    exi     = nan(2,D);
    for d=1:D
        mxi(d) = normx{d}(ahbi(d));
        exi(1,d) = normx{d}(ahbi(d)-sehbi(d));
        exi(2,d) = normx{d}(ahbi(d)+sehbi(d));
    end

    pxp{i}  = pxpi;
    Nbar{i} = Nbari;
    mnames{i} = mnamesi;
    pnames{i} = pnamesi;
    mx{i} = mxi;
    ex{i} = exi;
    
end

end