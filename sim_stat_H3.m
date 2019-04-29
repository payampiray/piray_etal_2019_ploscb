function [tr,xb,xl,xh,xc,ratehbi,ratebarhbi,ratesim,ratebarsim] = sim_stat_H3(fsimfits,m)

kref = 1;
ip   = 3;

nr  = length(fsimfits);
xb  = cell(nr,1);
xl  = cell(nr,1);
xh  = cell(nr,1);
xc  = cell(nr,1);

ratehbi = nan(nr,1);
% rate05lap = nan(nr,1);
ratesim = nan(nr,1);

ratebarhbi = nan(nr,2);
ratebarsim = nan(nr,2);

thr    = 0.05;

for i=1:nr    
    fsimfit = fsimfits{i};
%     [ahbi,sehbi,Nbar,xmsim,sesim,dfsim]=sim_stats_student(fsimfit,kref,ip);
    [ahbi,sehbi,Nbar,xmsim,sesim,dfsim]=sim_stats_student(fsimfit,kref,ip);
    
    nsamples = length(ahbi);
    

    thbi     = (ahbi- m )./sehbi;
    df       = Nbar+1;        
    pthbi    = 2 * tcdf(-abs(thbi), df);
%         rate05hbi(i)= mean(pthbi<thr);                
    [ratehbi(i),ratebarhbi(:,i)]=binofit(sum(pthbi<thr),nsamples);        

%     tlap     = (xmlap- m )./selap;
%     ptlap    = 2 * tcdf(-abs(tlap), dflap);
%     rate05lap(i)= mean(ptlap<thr);

    tsim     = (xmsim- m )./sesim;
    ptsim    = 2 * tcdf(-abs(tsim), dfsim);
%     ratesim(i)= mean(ptsim<thr);        
    [ratesim(i),ratebarsim(:,i)]=binofit(sum(ptsim<thr),nsamples); 
    
    m        = 0;
    t        = (ahbi- m )./sehbi;
    df       = Nbar+1;
    phbi     = tcdf(t, df);    
    [tr, xb{i},xl{i},xh{i},xc{i}] = sim_stats_student_bino(phbi,'cdf');
end

xb = cell2mat(xb);
xl = cell2mat(xl);
xh = cell2mat(xh);
xc = cell2mat(xc);

end
