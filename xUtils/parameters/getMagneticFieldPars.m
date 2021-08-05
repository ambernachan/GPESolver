function [p,q] = getMagneticFieldPars(Bz, Wmin, Ehfs)
    
    p = - getphysconst('ubohr') * Bz / (2*getphysconst('hbar')*Wmin);
    q = [];
    
    if nargin > 2
        q = getphysconst('ubohr')^2 * Bz^2 / (4*Ehfs*getphysconst('hbar')*Wmin);
    end
        
    
end