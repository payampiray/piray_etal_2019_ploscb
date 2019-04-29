function pnames = sim_plot_adjust_pnames(pnames)
% adjust pnames to be compatible with latex format

for k=1:length(pnames)
    if iscell(pnames{k})
        for i=1:length(pnames{k})
            pnames{k}{i} = tex2latex(pnames{k}{i});
        end
    else
        pnames{k} = tex2latex(pnames{k});
    end
end

end

function pname = tex2latex(pname)
if length(pname)>=4
    if strcmp(pname(1:4),'\it ')
        pname(1:4)='';
    end
end
pname = sprintf('$%s$',pname);
end