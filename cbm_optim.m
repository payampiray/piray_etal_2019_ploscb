function [tx,F,H,G,flag,k,P,NLL] = cbm_optim(h,pconfig,rng,numrep,init0,fid)
if nargin<4, numrep=1; end
if nargin<5, init0=[]; end
if nargin<6, fid  =1;  end

% P and NLL are vectors of prms and nll
options = optimset('LargeScale',pconfig.largescale,...
    'Display','off','TolFun',10^-10,'GradObj',...
    pconfig.gradient,'Hessian',pconfig.hessian);

rejcheck = 0;
if ~isempty(pconfig.reject_Fc)
    rejcheck = 1;
    srFc     = pconfig.reject_Fc;    
    rejectFc = str2func(srFc);
end
tolG       = pconfig.tolgrad;
numrep_up  = pconfig.numinit_up;
numrep_med = pconfig.numinit_med;
verbose    = pconfig.verbose;

F    = 10^16;
flag = 0;
tx   = nan;
H    = nan;
G    = nan;
r = rng(2,:)-rng(1,:); % for random initialization
numrep = numrep + size(init0,1);

if isempty(numrep_up), numrep_up = 100+numrep; end
k  = 0;
P  = [];
NLL= [];
while( (k<numrep) || (k>=numrep && k<numrep_med && flag==.5) || (k>=numrep && k<numrep_up && flag==0) )
    k=k+1;
    try
        init = init0(k,:);
    catch %#ok<CTCH>
        init = rand(size(r)).*r+rng(1,:);
    end
    try
        [tx_tmp, F_tmp, ~,output,grad_tmp,hess_tmp] = fminunc(h, init, options);
        [~,ishesspostmp] = chol(hess_tmp);
        ishesspostmp = ~logical(ishesspostmp);
        
        if rejcheck
            rejflag = rejectFc(F_tmp);
        else
            rejflag = 0;
        end
        sumG = mean(abs(grad_tmp));
        
        if (flag~=1 || (F_tmp<F)) && ishesspostmp && ~rejflag && (sumG<tolG)
                flag = 1;
                tx = tx_tmp;
                F = F_tmp;
                H = hess_tmp;
                G = grad_tmp;
        end
        if (flag~=1 && (F_tmp<F)) && ishesspostmp && ~rejflag && (sumG>tolG) % minimal condition
            flag  = .5;
            tx = tx_tmp;
            F  = F_tmp;
            H  = hess_tmp;
            G  = grad_tmp;
        end        
        
        P   = [P; tx_tmp]; %#ok<AGROW>
        NLL = [NLL; F_tmp]; %#ok<AGROW>
    catch msg
        fprintf(fid,'This initialization was aborted (there might be a problem with hfunc)\n');
        fprintf(fid,'The message of optimization routine is:\n') %#ok<PRTCAL>
        fprintf(fid,'   %s\n',msg.message);
    end
end

if verbose>0
    switch flag
        case 0
            fprintf(fid,'No positive hessian found in spite of %d initialization.\n',k);
        case .5
            fprintf(fid,'Positive hessian found, but not a good gradient in spite of %d initialization.\n',k);
            fprintf(fid,'We continue but be careful with this subject.\n');
        case 1
            if k>numrep
                fprintf(fid,'Optimized with %d initializations(>%d specified by user).\n',k,numrep);
            end
    end
end

end
