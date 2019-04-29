function [dx_hbi,dx_hier,dx_lap, mdxhbi, mdxhier, mdxlap] = sim_statfx_error(fit,normx)
dofx = 1;
if nargin<2, dofx = 0; fx=@(x)x; end

% xfit = fit.hbi.parameters;
lap  = fit.lap;
hier = fit.hier;
sz   = fit.sim.z;
xsim = fit.sim.h;
K    = length(lap);

[~,iz] = max(sz,[],1);

xfithbi  = cell(1,K);
xfitlap  = cell(1,K);
xfithier = cell(1,K);

for k=1:K
    xfithbi{k} = fit.hbi.parameters{k}';
    xfitlap{k} = lap(k).parameters';
    xfithier{k} = hier(k).parameters';
            
    np = size(lap(k).parameters,2);
    for i=1:np
        
        if dofx
            fx = normx{k}{i};
        end
        
        xfithbi{k}(i,:)=fx(xfithbi{k}(i,:));
        xfitlap{k}(i,:)=fx(xfitlap{k}(i,:));
        xfithier{k}(i,:)=fx(xfithier{k}(i,:));
        xsim{k}(i,:)=fx(xsim{k}(i,:));
    end
    
    
end

[dx_hbi , mdxhbi]  = statscomp(xfithbi,xsim,iz);
[dx_hier, mdxhier] = statscomp(xfithier,xsim,iz);
[dx_lap , mdxlap]  = statscomp(xfitlap,xsim,iz);

end

function [dx,mdx] = statscomp(xfit,xsim,im)
K = length(xsim);
dx = cell(1,K);

sdx   = 0;
nallp = 0;
for m=1:K
    imm  = im==m;    
    y    = abs(xfit{m}(:,imm)-xsim{m}(:,imm));
    dx{m}  = mean(y,2)';
    
    sdx    = sdx + sum(sum(y));
    nallp  = nallp+size(y,1)*sum(imm);
end
mdx = sdx/nallp;

end