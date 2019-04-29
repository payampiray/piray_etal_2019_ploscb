function [cmnew] = hbmc_cm(models,data,pconfig,qmutau,cm)
[N] = length(data);
[K] = length(models);

mu  = cell(K,1);
precision = cell(K,1);
priors = struct('mean',mu,'precision',precision);
D = nan(K,1);
for k=1:K
    D(k)      = length(qmutau(k).a);
    priors(k) = struct('mean',qmutau(k).a,'precision',diag(qmutau(k).Etau));
    for n=1:N
        inits = cm.theta{k}(:,n)';
        allconfigs(k,n) = pconfig(k); %#ok<AGROW>
        allconfigs(k,n).inits = inits; %#ok<AGROW>
    end
end

alg = pconfig(1).algorithm;

if strcmp(alg,'loophierlap')
    try
    [logf,theta,A,~,flag] = cbm_loop_quad_lap(data,models,priors,allconfigs,'hier');
    catch
        save('log_cm.mat','logf','theta','A','flag');
        error('error: %s\n',msg.message);        
    end
end

if strcmp(alg,'hierlap')
    theta = cell(K,1);
    logf  = nan(K,N);
    A     = cell(K,N);
    flag  = nan(K,N);
    for k=1:K
        for n=1:N
            [logf(k,n),theta{k}(:,n),A{k,n},~,flag(k,n)] = cbm_quad_lap(data{n},models{k},priors(k),pconfig(k),'hier'); 
        end
    end
end

Ainvdiag = cell(K,1);
logdetA  = nan(K,N);
for n=1:N  
    for k=1:K
        if flag(k,n)==0
            theta_n      = priors(k).mean;
            [logf(k,n)]  = cbm_loggaussian(theta_n',models{k},priors(k),data{n});    
            A{k,n}       = priors(k).precision;
            theta{k}(:,n)= theta_n;
        end
        logdetA_kn       = 2*sum(log(diag(chol(A{k,n}))));
        Ainvdiag{k}(:,n) = diag(A{k,n}^-1);
        logdetA(k,n)     = logdetA_kn;        
    end
end
cmnew   = struct('logf',logf,'theta',{theta},'Ainvdiag',{Ainvdiag},'logdetA',logdetA);
end