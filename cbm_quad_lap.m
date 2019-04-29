function [loglik,m,A,G,flag,kopt] = cbm_quad_lap(data,model,prior,config,mode)
% quadratic approximations using laplace fomula

switch mode
    case {'MAP','hier'}
        [loglik,m,A,G,flag,kopt] = quad_map(data,model,prior,config);
    case 'gaussian'
        [loglik,m,A,G,flag,kopt] = quad_gaussian(data,model,prior,config);
    otherwise
        error('mode %s is not recognized!',mode);
end
end

function [loglik,x,A,G,flag,kopt] = quad_map(data,model,prior,config)
pconfig = config;
pconfig.chance_v = prior.precision^-1;
% pconfig.verbose = 1;
numinit  = pconfig.numinit;
init0    = pconfig.inits;
rng      = pconfig.rng;

Gconf =  pconfig.gradient;
Hconf =  pconfig.hessian;

init  = [prior.mean'; init0];
hfunc = @(theta)(neglogloggaussian(theta,model,prior,data,Gconf,Hconf));
[tx,negloglik,H,G,flag,kopt] = cbm_optim(hfunc,pconfig,rng,numinit,init,1);
if flag == 0
    d  = length(prior.mean);
    tx = nan(1,d);
    H  = nan(d,d);
    G  = nan(d,1);
end
loglik = -negloglik;
x      = tx';
A      = H;
G      = G';
end

function [F,G,H] = neglogloggaussian(theta,model,prior,data,Gconf,Hconf)
[F,G,H] = cbm_loggaussian(theta,model,prior,data,Gconf,Hconf);
% % fun = @(theta)cbm_loggaussian(theta,model,prior,data,Gconf,Hconf);
% % [H2] = hessian(fun,theta);
% % [G2] = gradest(fun,theta);
F = -F;
G = -G;
H = -H;
end

function [loglik,tx,A,G,flag,kopt] = quad_gaussian(data,model,prior,config)

T  = prior.precision;
mu = prior.mean;

[~,PHI,t] = model(mu,data);

%
beta = 1;
Sinv = T + beta*(PHI'*PHI);
S    = Sinv^-1;
m    = S*(T*mu + beta*PHI'*t);

loglik = model(m,data);
G      = zeros(1,length(mu));
tx     = m;
A      = Sinv;
flag   = 1;
kopt   = 0;
end
