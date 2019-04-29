function math = hbmc_errorbar(math)
qmutau = math.qmutau;
K = length(qmutau);
for k=1:K
    nu    = qmutau(k).nu;
    beta  = qmutau(k).beta;
    sigma = qmutau(k).sigma;
    s2    = 2*sigma/beta;
    nk    = 2*nu;
    qmutau(k).nk = nk;
    qmutau(k).se = sqrt(s2/nk);
end
math.qmutau = qmutau;
end