function [value] = getphysconst(physconst)
    global hbar amu abohr kB c e me ubohr

    hbar = 1.0545718 * 10^(-34); % hbar in m^2 kg /s
    amu = 1.66053904 * 10^(-27); % Atom mass unit in kg
    abohr = 5.291772109 * 10^(-11); % Bohr radius in m
    kB = 1.38064852 * 10^(-23) ; % Boltzmann constant in kg m^2 / s^2 K
    c = 299792458; % speed of light (m/s)
    e = 1.60217662 * 10^(-19); % elementary charge in C (A*s)
    me = 9.1093897 * 10^(-31); % electron rest mass in kg
    ubohr = e*hbar / (2*me); % Bohr magneton in J/T or A*m^2
    
    if strcmp(physconst,'hbar')
        value = hbar;
    elseif strcmp(physconst,'amu')
        value = amu;
    elseif strcmp(physconst,'abohr')
        value = abohr;
    elseif strcmp(physconst,'kB')
        value = kB;
    elseif (strcmp(physconst,'c') || strcmp(physconst,'speedoflight'))
        value = c;
    elseif (strcmp(physconst,'e') || strcmp(physconst,'elementarycharge'))
        value = e;
    elseif (strcmp(physconst,'me') || strcmp(physconst,'electronmass'))
        value = me;
    elseif strcmp(physconst,'ubohr')
        value = ubohr;
    else
        error('Required physical constant is not recognized.')
    end
end