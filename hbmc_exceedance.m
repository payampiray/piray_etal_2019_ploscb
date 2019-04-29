function [hbmc,xp,pxp,bor] = hbmc_exceedance(hbmc,hbmc0)
if nargin<2, hbmc0 = []; end
inputisfile = 0;
if ischar(hbmc), inputisfile = 1; fhbmc=hbmc; hbmc = load(hbmc); hbmc = hbmc.math(end); end
if ischar(hbmc0), hbmc0 = load(hbmc0); hbmc0 = hbmc0.math(end); end

alpha = hbmc.qm.alpha;
L     = hbmc.bound.bound.L;
L0    = nan;
if isempty(hbmc0)
    if nargout>2
        error('For computing pxp and bor, hbmc0 is required!');
    end
else
    L0    = hbmc0.bound.bound.L;
end
Nsamp = 1e6;
[xp,pxp,bor] = compute_exceedance(alpha,L,L0,Nsamp);

hbmc.exceedance = struct('xp',xp,'pxp',pxp,'bor',bor,'alpha',alpha,'L',L,'L0',L0,'hbmc0',hbmc0);
if inputisfile, save(fhbmc,'hbmc'); end
end

function [xp,pxp,bor] = compute_exceedance(alpha,L,L0,Nsamp)
K     = length(alpha);
if K == 2
    % comparison of 2 models
    xp(1) = betacdf(0.5,alpha(2),alpha(1));
    xp(2) = betacdf(0.5,alpha(1),alpha(2));
else
    % comparison of >2 models: use sampling approach
    xp = dirichlet_exceedance(alpha,Nsamp);
end

bor = 1/(1+exp((L-L0)));
% Compute protected exceedance probs - Eq 7 in Rigoux et al.
pxp=(1-bor)*xp+bor/K;

end


function xp = dirichlet_exceedance(alpha,Nsamp)
% (almost) entirely copied from spm toolbox, 2008 Wellcome Trust Centre for Neuroimaging
% spm_dirichlet_exceedance
% version:
% $Id: spm_dirichlet_exceedance.m 3118 2009-05-12 17:37:32Z guillaume $
% see also spm_BMS
 
% Compute exceedance probabilities for a Dirichlet distribution
% FORMAT xp = spm_dirichlet_exceedance(alpha,Nsamp)
% 
% Input:
% alpha     - Dirichlet parameters
% Nsamp     - number of samples used to compute xp [default = 1e6]
% 
% Output:
% xp        - exceedance probability
%__________________________________________________________________________
%
% This function computes exceedance probabilities, i.e. for any given model
% k1, the probability that it is more likely than any other model k2.  
% More formally, for k1=1..Nk and for all k2~=k1, it returns p(x_k1>x_k2) 
% given that p(x)=dirichlet(alpha).
% 
% Refs:
% Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ
% Bayesian Model Selection for Group Studies. NeuroImage (in press)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny & Klaas Enno Stephan
% $Id: spm_dirichlet_exceedance.m 3118 2009-05-12 17:37:32Z guillaume $

if nargin < 2
    Nsamp = 1e6;
end

Nk = length(alpha);

% Perform sampling in blocks
%--------------------------------------------------------------------------
blk = ceil(Nsamp*Nk*8 / 2^28);
blk = floor(Nsamp/blk * ones(1,blk));
blk(end) = Nsamp - sum(blk(1:end-1));

xp = zeros(1,Nk);
for i=1:length(blk)
    
    % Sample from univariate gamma densities then normalise
    % (see Dirichlet entry in Wikipedia or Ferguson (1973) Ann. Stat. 1,
    % 209-230)
    %----------------------------------------------------------------------
    r = zeros(blk(i),Nk);
    for k = 1:Nk
%         r(:,k) = spm_gamrnd(alpha(k),1,blk(i),1);
        r(:,k) = gamrnd(alpha(k),1,blk(i),1);
%             r = gamrnd(repmat(vec(suffStat.d),1,N),1,K,N);
%             y = r ./ repmat(sum(r,1),K,1);        
    end
    sr = sum(r,2);
    for k = 1:Nk
        r(:,k) = r(:,k)./sr;
    end
    
    % Exceedance probabilities:
    % For any given model k1, compute the probability that it is more
    % likely than any other model k2~=k1
    %----------------------------------------------------------------------
    [~, j] = max(r,[],2);
    xp = xp + histc(j, 1:Nk)';
    
end
xp = xp / Nsamp;

end
