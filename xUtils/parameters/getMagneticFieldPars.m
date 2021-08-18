function [p,q] = getMagneticFieldPars(Bz, Wmin, Ehfs)
    
    p = - getphysconst('ubohr') * Bz / (2*getphysconst('hbar')*Wmin);
    % From Stamper-Kurn, Ueda I get a 100*2pi difference in the linear
    % Zeeman energy parameter (compared to Bao's paper), here's the diff:
    p = p / (2*pi*100);
    q = [];
    
    if nargin > 2
        q = getphysconst('ubohr')^2 * Bz^2 / (4*Ehfs*getphysconst('hbar')*Wmin);
    end
        
    
end