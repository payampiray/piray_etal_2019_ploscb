function [bound,dL] = hbmc_bound(bound,lastmodule)
bb            = bound.bound;
pmlimInf      = bb.pmlimInf;
Elogpm_Elogqm = bb.Elogpm - bb.Elogqm;
if pmlimInf
    Elogpm_Elogqm = 0;
end
Lpre     = + bb.ElogpX  + bb.ElogpH   + bb.ElogpZ + ...
           + bb.Elogpmu + bb.Elogptau + ...
           - bb.ElogqH  - bb.ElogqZ +...
           - bb.Elogqmu - bb.Elogqtau + ...
           + Elogpm_Elogqm;
       
switch lastmodule
    case 'qHZ'
        bb.ElogpX    = sum(bound.qHZ.ElogpX);
        bb.ElogpH    = sum(bound.qHZ.ElogpH);
        bb.ElogpZ    = sum(bound.qHZ.ElogpZ);
        bb.ElogqH    = sum(bound.qHZ.ElogqH);
        bb.ElogqZ    = sum(bound.qHZ.ElogqZ);        
        
    case 'qmutau'
        
        bb.ElogpH    = sum(bound.qmutau.ElogpH);
        bb.Elogpmu   = sum(bound.qmutau.Elogpmu);
        bb.Elogptau  = sum(bound.qmutau.Elogptau);
        bb.Elogqmu   = sum(bound.qmutau.Elogqmu);
        bb.Elogqtau  = sum(bound.qmutau.Elogqtau);
        
    case 'qm'
        bb.ElogpZ    = sum(bound.qm.ElogpZ);
        bb.Elogpm    = bound.qm.Elogpm;
        bb.Elogqm    = bound.qm.Elogqm;
end
Elogpm_Elogqm = bb.Elogpm - bb.Elogqm;
if pmlimInf
    Elogpm_Elogqm = 0;
end
L        = + bb.ElogpX  + bb.ElogpH   + bb.ElogpZ + ...
           + bb.Elogpmu + bb.Elogptau + ...
           - bb.ElogqH  - bb.ElogqZ +...
           - bb.Elogqmu - bb.Elogqtau + ...
           + Elogpm_Elogqm;
dL       = + L - Lpre;


bb.lastmodule = lastmodule;
bb.L          = L;
bb.dL         = dL;
bound.bound   = bb;
bound.(lastmodule).bound = bb;

end