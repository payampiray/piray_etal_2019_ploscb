function cbm = cbm_lap_aggregate(cbms)

N         = length(cbms);
data      = cell(N,1);

flags     = nan(1,N);
telapsed  = nan(1,N);

loglik    = nan(1,N);
theta     = cell(1,N);
T         = cell(1,N);
lme       = nan(1,N); % log-model-evidence
Tinv_diag = cell(1,N);

for n=1:N
    cbm = load(cbms{n}); cbm = cbm.cbm;
    % input field
    if ~isempty(cbm.input.data)
        data(n)  = cbm.input.data;
    end
    if n==1
        model          = cbm.input.model;
        functionname   = cbm.input.functionname;
        prior          = cbm.input.prior;
        config         = cbm.input.config;
        fname          = [];
    end
    % profile field
    flags(n)     = cbm.profile.optim.flag;
    gradient(:,n)= cbm.profile.optim.gradient; %#ok<AGROW>
    telapsed(n)  = cbm.profile.optim.telapsed;
    if n==1
        numinit  = cbm.profile.optim.numinit;
        rng      = cbm.profile.optim.rng;
    end
    
    % math field
    loglik(n)    = cbm.math.loglik;
    theta(n)     = cbm.math.theta;
    T(n)         = cbm.math.T;
    lme(n)       = cbm.math.lme;
    Tinv_diag(n) = cbm.math.Tinv_diag;
    
    % output field
    lme(n) = cbm.output.log_evidence;
end

input      = struct('data',{data},'model',model,...
                    'functionname',functionname,'prior',prior,'config',config,'fname',fname);
math       = struct('loglik',loglik,'theta',{theta},'T',{T},'lme',lme,'Tinv_diag',{Tinv_diag});
optim      = struct('numinit',numinit,'rng',rng,'telapsed',telapsed,'flag',flags,'gradient',gradient);
profile    = struct('datetime',datestr(now),'filename',mfilename,'optim',optim);
output     = struct('parameters',cell2mat(theta)','loglik',loglik','log_evidence',lme');
cbm        = struct('method',mfilename,...
                    'input',input,...
                    'profile',profile,...
                    'math',math,...                    
                    'output',output);
end