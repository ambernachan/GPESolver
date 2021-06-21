function [XI] = findhealinglengths(chis, scatlengths, atom)

    % required parameters / constants
    N = getsimconst('N'); % number of particles
    hbar = getphysconst('hbar'); % in kg m^2 / s
    trapfreq = getsimconst('trap_freq'); % Trap strength in Hz (symmetric for now)
    atom_mass = getsimconst(['mass_' atom]); % Atom mass in kg
    aho = sqrt(hbar / (atom_mass * trapfreq)); % characteristic length osc in m

    XI = zeros(1,length(chis));
    for i=1:length(chis)
        chi = chis(i);
        a = scatlengths(i);
        if abs(chi) >= 10 % strong interactions
            xi = abs(aho * ( aho / (15*N*a) )^(1/5));
        elseif abs(chi) <= 0.1 && ~(abs(chi) < 0.01) % weak interactions !=0
            xi = 1 / sqrt(8*pi*a*getsimconst('density'));
        elseif abs(chi) == 0 || (abs(chi) < 0.01) % (almost) no interactions
            xi = aho;
        elseif abs(chi) > 0.1 && abs(chi) < 10 % in-between state
            xi = abs(0.5*aho * ( aho / (15*N*a) )^(1/5)) ...
                + abs(0.5 / sqrt(8*pi*a*getsimconst('density')));
        end
        XI(i) = xi;
    end
    str = ['Healing length for self-interactions: xi ~ %0.0f nm\n' ...
          'Healing length for spin-mixing interactions: xi ~ %0.0f nm'];
    sprintf(str, XI(1)*10^9, XI(2)*10^9)
end