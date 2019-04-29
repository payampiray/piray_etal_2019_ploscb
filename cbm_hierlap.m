function cbm = cbm_hierlap(data, model, fcbm_map, fname, pconfig)
% hierarchical-laplace approximation
% [CBM, SUCCESS] = cbm_hierlap(DATA, MODEL, FCBM_MAP, FNAME, PCONFIG)
% DATA:  data (Nx1 cell) where N is number of samples
% MODEL: a function-handle to the model computing log-likelihood of DATA 
% given some paramete 
% FCBM_MAP: file-address of a maximum-a-posteriori approximation for each
% subject in the cbm-compatible format. This fle could be for example the
% output of cbm_lap.
% FNAME: filename for saving the output (leave it empty for not saving)
% PCONFIG: a struct for configuration (leave it empty if you do not have
% any idea!), see cbm_config for more info
% CBM: the main output-file with following fields:
%   CBM.output is output field. 
%       CBM.output.parameters is N-by-d matrix containing fitted parameters.
%       CBM.output.group_mean is 1-by-d containing mean of group
%       parameters.
%       CBM.output.group_variance is 1-by-d containing variance of group
%       parameters.
%       CBM.output.log_evidence is a scaler representing log-model evidence. 
%   CBM.input contains inputs to cbm_hierlap
%   CBM.profile contains info about optimization
%   CBM.math contains all relevant variables 
% see cbm_example_hierlap for an example.
% 
% This function performs a hierarchical model-fitting using laplace 
% approximation of individual posteriors. This function needs that every 
% subject had been fitted independently before (for example by cbm_lap).
% Model-fitting is based on a expectation-maximization procedure, in which
% individual parameters are optimized in the (laplace approximated) E-step 
% and group mean and variance are optimized in the M-step 
% Model-evidence is based on Laplace approximation for individual
% parameters and Bayesian-information-criterion for 2d group-level
% parameters
% 
% Important configurations
% tolx is tolerance for changes in (normalized valeu of) mean parameters; 
% defaults is 0.0125=norminv(0.505)
% tolL is tolerance for changes in the lower bound; default is 0.05
% maxiter indicates the maximum number of iterations; default is 50
% terminate (dx or dL) indicates termination variable. dx indicates that
% the algorithm stops if dx is lower than tolx. dL indicates that the
% algorithm stops if dL is lower than tolL.
% 
% Other configurations
% verbose: whether to write on screen (or a log-file), default is 1
% functionname: an arbitray name (for example the model-name), default is
% empty
% fname_prog: a file-address for saving progress (iterations) is save_prog
% is one, default cbm_[NOW].mat in the current directory
% rng: a 2-by-d matrix representing the range of parameters. random seeds
% for initializations will be picked in this range; default 
% [-5*ones(1,d);5*ones(1,d)]
% save_data: a 0 or 1 scaler indicating whether to save data in the output
% file or not; default is 1 (save data)
% numinit: number of random initialization, default is 0
% save_prog: a 0 or 1 scaler indicating whether to save progress (in the
% file indicated by fname_prog) or not; default is 1.
% numinit_med is the (increased) number of random initializations for bad
% subjects; default is 10.
% numinit_med is the (maximum) number of random initializations for bad
% subjects; default is 50.
% 
% dependecies:
% this function needs MATLAB optimization toolbox (calls fminunc)
% cbm_config, cbm_check_input, cbm_quad_lap, cbm_optim, cbm_loop_quad_lap
% 
% cite this paper if you use this function for model-comparison: 
% Piray et al., 2014, J Neurosci
% cite paper if you use this function for model-fitting:
% Huys et al., 2011, PLoS Comp Biol
% 
% Note: this code is free (under GNU open-source licence), but you use 
% it with your own responsibility!
% [if you found any bug, please let me know]
% April 2017, Payam Piray 
%==========================================================================

%----------------------------------
% initialize using cbm_map
fcbmerrmsg = 'fcbm_map input has not properly been specified!';
if ischar(fcbm_map) % address of a cbm saved by cbm_lap?
    cbm_map  = load(fcbm_map); cbm_map = cbm_map.cbm;
elseif isstruct(fcbm_map) % or itself a cbm struct?
    cbm_map  = fcbm_map;
else % or something is wrong
    error(fcbmerrmsg);
end
try
    theta0   = cbm_map.math.theta;
    logf0    = cbm_map.math.loglik;
    T        = cbm_map.math.T;
    Tinvdiag = arrayfun(@(n)diag(T{n}^-1),(1:length(T)),'UniformOutput',0);
    logdetT  = cell2mat(arrayfun(@(n)log(det(T{n})),(1:length(T)),'UniformOutput',0));
    Gr       = arrayfun(@(n)cbm_map.profile.optim.gradient(:,n)',(1:length(T)),'UniformOutput',0);
    flag     = cbm_map.profile.optim.flag;
    prior0   = cbm_map.input.prior;
catch
    error(fcbmerrmsg);
end
quad_apx = struct('logf',logf0,'theta',{theta0},'T',{T},'Tinvdiag',{Tinvdiag},'logdetT',logdetT);
d        = size(theta0{1},1);

T0 = prior0.precision;
if ~isvector(T0)
    V0 = diag(T0^-1);
else
    V0 = T0.^-1;
end
prior0.variance = V0;

%----------------------------------
% fitting specs
if nargin<4, error('Not enough inputs'); end
if nargin<5, pconfig=[]; end

looplap     = 0;
if ~isfield(pconfig,'algorithm')
    algorithm = 'hierlap';
else
    algorithm = pconfig.algorithm;
end
if strcmp(algorithm,'loophierlap')
    looplap = 1;
end

pconfig     = cbm_config(d,algorithm,pconfig);
fnameprog   = pconfig.fname_prog;
saveprog    = pconfig.save_prog;
numinit     = pconfig.numinit;
flog        = pconfig.flog;
verbose     = pconfig.verbose;


freemu      = pconfig.free_group==1;
freevar     = pconfig.free_groupvar==1;

% functionname, for display
functionname = pconfig.functionname;
if isempty(functionname), functionname = func2str(model); end

%----------------------------------
% log?
fid  = 1;
if ischar(flog)
    fid = fopen(flog,'w');
end

%----------------------------------
% check inputs
ok = cbm_check_input(rand(d,1)',model,data,fname);
if ok
    if verbose, fprintf(fid,'Model is good! Continue with fitting...\n'); end
else
    if fid~=1, fprintf(fid,'Model is not good at mean prior! Have to stop here, sorry!\n'); end
    error('Model is not good at mean prior! Have to stop here, sorry!\n');
end
if verbose, fprintf(fid,'%-70s\n',repmat('-',1,70)); end

%----------------------------------
% intro!
tic;
if verbose, fprintf(fid,'%-40s%30s\n',mfilename,datestr(now)); end
if verbose, fprintf(fid,'%-70s\n',repmat('=',1,70)); end

%----------------------------------
% dimensions
N   = length(data); % number of samples (subjects)
if verbose, fprintf(fid,'function to be optimized: %s\n',functionname); end
if verbose, fprintf(fid,'Number of samples: %d\n',N); end
if verbose, fprintf(fid,'Number of parameters: %d\n\n',d); end
if verbose, fprintf(fid,'Number of initializations: %d\n',numinit+2); end
if verbose, fprintf(fid,'%-70s\n',repmat('-',1,70)); end

%==========================================================================
        
%----------------------------------
% iterations config
dx        = nan;
dL        = nan;
L_pre     = nan;
mu_pre    = nan(d,1);
tolx      = pconfig.tolx;
tolL      = pconfig.tolL;
iter      = 0;
maxiter   = pconfig.maxiter;

mu = 0;
W  = 0;
%----------------------------------
% termination
prog     = struct('L',[],'dL',[],'dx',[]);
dprog    = struct('dx',dx,'dL',dL,'L_pre',L_pre,'mu_pre',mu_pre);
switch pconfig.terminate
    case 'dx'
        terminate = @(dx,dL)(dx<tolx);
    case 'dL'
        terminate = @(dx,dL)(abs(dL)<tolL);
    case 'dxdL'
        terminate = @(dx,dL)(abs(dL)<tolL) && (dx<tolx);
end

%----------------------------------
% main algorithm
while  ~terminate(dprog.dx,dprog.dL) && (iter<maxiter)
    iter       = iter + 1;
    
    if iter>1
        logf       = nan(1,N);
        logdetT    = nan(1,N);
        theta      = theta0;
        T          = cell(1,N);
        Tinvdiag  = cell(1,N);
        Gr         = cell(1,N);
        flag       = nan(1,N);    
        prior  = struct('mean',mu,'precision',diag(W.^-1));    
        if verbose, fprintf(fid, 'Iteration %02d',iter); end
        if ~looplap
            if verbose, fprintf(fid, ', sample '); end
            for n=1:N        
                if verbose, fprintf(fid, '%02d/%02d',n,N); end
                % -data
                dat = data{n}; 

                % -initialization: (mu and) the previous one        
                inits(1,:) = theta{n}';

                % -optimize
                pconfig.inits = inits;
                [logf_n,theta_n,T_n,Gr_n,flag_n] = cbm_quad_lap(dat,model,prior,pconfig,'hier');        

                % -check if optimization ok, otherwise use prior values
                if flag_n==0
                    if verbose, fprintf(fid,'oops! bad subject... use group parameters for this one.\n'); end
                    theta_n   = mu;
                    [logf_n]  = cbm_loggaussian(theta_n',model,prior,dat);
                    Gr_n      = nan(1,d);
                    T_n       = T0;
                end

                logf(n)       = logf_n;
                theta{n}      = theta_n;
                T{n}          = T_n;
                Gr{n}         = Gr_n;
                Tinvdiag{n}   = diag(T_n^-1);
                logdetT(n)    = 2*sum(log(diag(chol(T_n)))); %=log(det(T_n));                
                flag(n)       = flag_n;
                if verbose && (n<N) && fid==1, fprintf(fid, '\b\b\b\b\b'); end
            end
        else
            pconfigloop(1:N) = deal(pconfig);
            for n=1:N
                inits(1,:) = theta{n}';
                pconfigloop(n).inits = inits;
            end
            [logf,thetaall,T,Gr,flag] = cbm_loop_quad_lap(data,{model},prior,pconfigloop,'hier');  
            for n=1:N   
                if flag(n)==0
                    thetaall{1}(:,n)  = priors(k).mean;
                    [logf(k,n)]       = cbm_loggaussian(thetaall{1}(:,n)',model,prior,data{n});    
                    T{n}              = T0;                    
                end
                theta{n}    = thetaall{1}(:,n);
                Tinvdiag{n} = diag(T{n}^-1);
                logdetT(n)  = 2*sum(log(diag(chol(T{n}))));
            end
        end
        quad_apx = struct('logf',logf,'theta',{theta},'T',{T},'Tinvdiag',{Tinvdiag},'logdetT',logdetT);
        if verbose, fprintf(fid, '...completed\n'); end    
    end
    
    % update group parameters
    sumstat        = cal_sumstat(quad_apx,freemu,freevar,prior0);    
    
    % get ready for the next iteration
    L              = sum(quad_apx.logf +d/2*log(2*pi) -quad_apx.logdetT/2);
    mu             = sumstat.mu;
    W              = sumstat.V;
    [prog, dprog]  = cal_converge(prog,dprog,L,mu,W);
    
    % wrap-up this iteration
    math        = struct('quad_apx',quad_apx,'sumstat',sumstat,'bound',L,'initialize',cbm_map.math);
    cbm_i(iter) = struct('math',math,'prog',prog,'dprog',dprog,'flag',flag,'G',{Gr});     %#ok<AGROW>
    if saveprog, save(fnameprog,'cbm_i'); end
       
    if verbose, fprintf(fid,'%-40s%30s\n',' ',sprintf('dx: %7.5f',dprog.dx)); end
    if verbose, fprintf(fid,'%-40s%30s\n',' ',sprintf('dL: %+7.3f',dprog.dL)); end
end

%----------------------------------
% group model evidence: penalty for 2d group free parameters
dg = sum(freemu) + sum(freevar);
group_lme  = math.bound -.5*dg*log(N);

%----------------------------------
% output
sdata      = {}; if pconfig.save_data, sdata = data; end
input      = struct('data',{sdata},'model',func2str(model),...
                    'functionname',functionname,'fname_cbm_map',fcbm_map,'config',pconfig,'fname',fname);

cbm_ib     = cbm_i(end);
math       = cbm_ib.math;
prog       = cbm_ib.prog;
theta      = cell2mat(math.quad_apx.theta);
mu         = math.sumstat.mu;
V          = math.sumstat.V;
output     = struct('parameters',theta','group_mean',mu','group_variance',V','log_evidence',group_lme,'free_group',freemu,'free_groupvar',freevar);
optim      = struct('numinit',numinit,'flag',cbm_ib.flag,'gradient',{cbm_ib.G});
profile    = struct('datetime',datestr(now),'filename',mfilename,'optim',optim);    

cbm        = struct('method','hierlap',...
                    'input',input,...
                    'profile',profile,...
                    'math',math,...
                    'prog',prog,...
                    'output',output);
if ~isempty(fname), save(fname,'cbm'); end

end

function sumstat = cal_sumstat(quad_apx,freemu,freevar,prior0)
mu0         = prior0.mean;
V0          = prior0.variance;

theta       = cell2mat(quad_apx.theta);
Tinvdiag    = cell2mat(quad_apx.Tinvdiag);
mu          = mean(theta,2);
mu(~freemu) = mu0(~freemu);

Vn          = bsxfun(@minus,theta,mu).^2;
V           = mean(Vn+Tinvdiag,2);

V(~freevar)  = V0(~freevar);

sumstat     = struct('mu',mu,'V',V);
end

function [prog, dprog]  = cal_converge(prog,dprog,L,mu,V)
L_pre  = dprog.L_pre;
mu_pre = dprog.mu_pre;

dL      = L - L_pre;
L_pre   = L;
mu_norm = mu./sqrt(V);
dx      = sqrt(mean((mu_norm-mu_pre).^2));

prog.L   = [prog.L L];
prog.dL  = [prog.dL dL];
prog.dx  = [prog.dx dx];
dprog    = struct('dx',dx,'dL',dL,'L_pre',L_pre,'mu_pre',mu_norm);
end

