function [value] = getphysconst(physconst)
    global hbar amu abohr kB c

    hbar = 1.0545718 * 10^(-34); % hbar in m^2 kg /s
    amu = 1.66053904 * 10^(-27); % Atom mass unit in kg
    abohr = 5.291772109 * 10^(-11); % Bohr radius in m
    kB = 1.38064852 * 10^(-23) ; % Boltzmann constant in kg m^2 / s^2 K
    c = 299792458; % speed of light (m/s)
    
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
    else
        error('Required physical constant is not recognized.')
    end
end