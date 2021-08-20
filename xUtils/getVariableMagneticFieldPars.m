function [P,Q] = getVariableMagneticFieldPars(Bmax, Bmin, Wmin, Ehfs, xlim)
    
    BZ = @(z) Bmin + (Bmax-Bmin)*(1+z/xlim)/2;
    
    % From Stamper-Kurn, Ueda I get a 100*2pi difference in the linear
    % Zeeman energy parameter (compared to Bao's paper), here's the diff:
    P = @(z) - getphysconst('ubohr') * BZ(z) / (2*getphysconst('hbar')*Wmin) / (2*pi*100);
    Q = @(z) 0;

    if nargin > 3
        Q = @(z) getphysconst('ubohr')^2 * BZ(z).^2 / (4*Ehfs*getphysconst('hbar')*Wmin);
    end
    
end