function [value] = getsimconst(simconst)
    % Returns simulations constants. Included are: a0 and a2 for Na and Rb;
    % mass_Na and Rb; number of particles N; trap frequency, spin pair no.

    a0_Na = 48.90 * getphysconst('abohr'); % scattering length (m)
    a2_Na = 54.54 * getphysconst('abohr'); % scattering length (m)
    a0_Rb = 101.65 * getphysconst('abohr'); % scattering length (m)
    a2_Rb = 100.49 * getphysconst('abohr'); % scattering length (m)
    
    mass_Na = 22.989769 * getphysconst('amu'); % mass (kg)
    mass_Rb = 86.9091835 * getphysconst('amu'); % mass (kg)
    rvdW_Na = 44.96 * getphysconst('abohr'); % van der Waals radius (m)
    rvdW_Rb = 82.64 * getphysconst('abohr'); % van der Waals radius (m)
    TNa = 3.24420293315 * 10^(-7); % ?
    TRb = 8.58177228994 * 10^(-8); % ?
    
    N = 10^5; % number of particles
    trap_freq = 100 * 2*pi; % in Hz
    spin_pair = 1; % hyperfine spin manifold
    
    ws_variables = whos;
    wsnames = cell(1,length(ws_variables));
    for i = 1:length(ws_variables)
        wsnames{i} = ws_variables(i).name;
    end
    
    for i = 1:length(wsnames)
        if strcmp(simconst, wsnames{i})
            value = eval(wsnames{i});
        end
    end
    
    % Check whether 'value' now exists
    if ~exist('value','var')
        error('Unrecognized simulation constant required.')
    end
        

end