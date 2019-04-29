function [tr,xb,xl,xh,xc] = sim_stat_H2(fsimfits,m)
[nr] = length(fsimfits);

kref = 1;
ip   = 3;


xb  = cell(nr,1);
xl  = cell(nr,1);
xh  = cell(nr,1);
xc  = cell(nr,1);

for i=1:nr    
    
    fsimfit = fsimfits{i};        
    [ahbi,sehbi,Nbar]=sim_stats_student(fsimfit,kref,ip);    
    
    thbi     = (ahbi- m )./sehbi;
    df       = Nbar+1;
    phbi     = tcdf(thbi, df);    
        
    [tr, xb{i},xl{i},xh{i},xc{i}] = inspect_binoci(phbi,'cdf');                
end

xb = cell2mat(xb);
xl = cell2mat(xl);
xh = cell2mat(xh);
xc = cell2mat(xc);

end

function [trlabels, xb,xl,xh,xc] = inspect_binoci(p,mode)
nsamples = length(p);

% [p,I] = sort(p,2);
% df = df(I);
% t  = t(I);

N = 20;
step = 1/N;

if strcmp(mode,'pdf')
tr = step:step:(1);
trpre = 0;
pdfb = zeros(1,N);
pdfl = zeros(1,N);
pdfh = zeros(1,N);
pdfd = zeros(1,N);
trlabels = cell(1,N);
for j=1:length(tr)
    pp = (p<=tr(j)) & (p>=trpre);   
    nj     = sum(pp ,2);%+ sum(tp{k}>0.95 ,2);
%     if numsim==2, nj = binornd(nsamples,0.05); end
    [bjp,bjci]=binofit(nj,nsamples);
    pdfb(1,j) = bjp; 
    pdfl(1,j) = bjci(1); 
    pdfh(1,j) = bjci(2); 
    pdfd(1,j) = diff(bjci);
%     pdfc(1,j) = ((step>=bjci(1)) && (step<=bjci(2))) + 0.000000001;
    trlabels{j} = sprintf('%0.2f-%0.2f',trpre,tr(j));
    
    trpre= tr(j);  
end

res  = 1/nsamples;
pdfl = floor(pdfl*res^-1)*res;
pdfh = ceil(pdfh*res^-1)*res;
pdfc = (pdfh>=step) & (pdfl<=step);

xb = pdfb;
xc = pdfc;
xl = pdfl;
xh = pdfh;

end

%-----------
% CDF

if strcmp(mode,'cdf')

tr = 0:step:1; trpre = 0;
N = length(tr);
cdfb = nan(1,N);
cdfl = nan(1,N);
cdfh = nan(1,N);
cdfd = nan(1,N);
for j=1:length(tr)
    pp = (p<=tr(j)) & (p>=trpre);
    nj     = sum(pp ,2);%+ sum(tp{k}>0.95 ,2);
%     if numsim==2, nj = binornd(nsamples,tr(j)); end
    [bjp,bjci]=binofit(nj,nsamples);
    cdfb(1,j) = bjp; 
    cdfl(1,j) = bjci(1); 
    cdfh(1,j) = bjci(2); 
    cdfd(1,j) = diff(bjci);
%     cdfc(1,j) = (tr(j)>bjci(1)) && (tr(j)<bjci(2));
end

res  = 1/nsamples;
cdfl = floor(cdfl*res^-1)*res;
cdfh = ceil(cdfh*res^-1)*res;
cdfc = (cdfh>=tr) & (cdfl<=tr);

xb = cdfb;
xc = cdfc;
xl = cdfl;
xh = cdfh;
trlabels = tr;
end

end
