function [rate05,methods] = sim_stat_H1(fsimfits,m)
[nr, nc] = size(fsimfits);

kref = 1;
ip   = 3;

thrhier = 3;
thr = .05;

rate05  = cell(nr,1);
for i=1:nr    
    rlap  = nan(1,nc);
    rhier = nan(1,nc);
    rhbi  = nan(1,nc);
    for j=1:nc
        fsimfit = fsimfits{i,j};        
        [ahbi,sehbi,Nbar,xmsim,sesim,dfsim,xmlap,selap,dflap,Fhier]=sim_stats_student(fsimfit,kref,ip);
        
        thbi     = (ahbi- m )./sehbi;
        df       = Nbar+1;        
        phbi     = 2 * tcdf(-abs(thbi), df);
        
        tlap     = (xmlap- m )./selap;
        plap     = 2 * tcdf(-abs(tlap), dflap);

        tsim     = (xmsim- m )./sesim;
        psim    = 2 * tcdf(-abs(tsim), dfsim);        

        
        rlap(j)  = corrclassify(psim<=thr,plap<=thr);
        rhier(j) = corrclassify(psim<=thr,Fhier>=thrhier);
        rhbi(j)  = corrclassify(psim<=thr,phbi<=thr);        
    end
    
    rate05{i}  = [rlap; rhier; rhbi];
end
methods = {'NHI','HPE','HBI'};


end

function [rate05,power,specificity] = corrclassify(p0,p)

TP = sum( p0 & p);
FN = sum( p0 & ~p);    
FP = sum( ~p0 & p);
TN = sum( ~p0 & ~p);

% TPR = TP/(TP+FN);
% TNR = TN/(TN+FP);
% 
% FPR = FP/(TN+FP);

% Power (equal to 1 - type II error)
beta = FN/(TP+FN);
power = 1-beta;

% Specificity or 1 - type I error
alpha = FP/(FP+TN);
specificity = alpha;


% TR  = (TP+TN)/(TP+TN+FP+FN);
rate05 = mean(p);
end