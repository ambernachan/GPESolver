function [value] = getsimconst(simconst)
    % Returns simulations constants. Included are: a0 and a2 for Na and Rb;
    % mass_Na and Rb; number of particles N; trap frequency, spin pair no.
    
    axicon_radius = 288.75 * 10^(-6); % axicon radius (m)
    laser_power = 1.5; % laser power (W)
    laser_freq = 563.53 * 2*pi * 10^(12); % laser frequency (Hz)
    laser_wavelength = 532 * 10^(-9); % laser wavelength (m)
    dipole_waist_x = 31.26 * 10^(-6); % dipole beam waist in x (m)
    dipole_waist_y = 406.42 * 10^(-6); % dipole beam waist in x (m)

    a0_Na = 48.90 * getphysconst('abohr'); % scattering length (m)
    a2_Na = 54.54 * getphysconst('abohr'); % scattering length (m)
    a0_Rb = 101.65 * getphysconst('abohr'); % scattering length (m)
    a2_Rb = 100.49 * getphysconst('abohr'); % scattering length (m)
    a0_ferro = 100 * getphysconst('abohr'); % scattering length (m)
    a2_ferro = 10 * getphysconst('abohr'); % scattering length (m)
    
    mass_Na = 22.989769 * getphysconst('amu'); % mass (kg)
    mass_Rb = 86.9091835 * getphysconst('amu'); % mass (kg)
    transitionfreq_Na = 508.8487162 * 2*pi * 10^(12); % transition frequency (Hz)
    linewidth_Na = 9.7946 *2*pi * 10^6; % natural linewith (Hz)
    transitionfreq_Rb = transitionfreq_Na;
    linewidth_Rb = linewidth_Na;
    Ehfs_Na = 1771.6261288 * 2*pi * 10^(6) * getphysconst('hbar'); % hyperfine energy splitting (J)
    Ehfs_Rb = 6835 * 2*pi * 10^(6) * getphysconst('hbar'); % hyperfine energy splitting (J)
    
    mass_ferro = 0.5*mass_Na + 0.5*mass_Rb;
    rvdW_Na = 44.96 * getphysconst('abohr'); % van der Waals radius (m)
    rvdW_Rb = 82.64 * getphysconst('abohr'); % van der Waals radius (m)
    TNa = 3.24420293315 * 10^(-7); % ?
    TRb = 8.58177228994 * 10^(-8); % ?
    
    N = 2 * 10^6; % number of particles
    trap_freq = 1 * 2*pi; % in Hz (minimum trap freq)
    spin_pair = 1; % hyperfine spin manifold
    density = 10^18; % density n in 1/m^3
    
    detuning_Na = - abs(laser_freq - transitionfreq_Na);
    detuning_Rb = - abs(laser_freq - transitionfreq_Rb);
    zRx = dipole_waist_x^2 * pi / laser_wavelength;
    zRy = dipole_waist_y^2 * pi / laser_wavelength;
    
    dipoleTrap0_unscaled_Na = 6*pi^2*getphysconst('c')^2 * laser_power * linewidth_Na ...
        / (transitionfreq_Na^2 * laser_wavelength * detuning_Na * (transitionfreq_Na + laser_freq)) ...
        / sqrt(zRx*zRy);
    dipoleTrap0_unscaled_Rb = 6*pi^2*getphysconst('c')^2 * laser_power * linewidth_Rb ...
        / (transitionfreq_Rb^2 * laser_wavelength * detuning_Rb * (transitionfreq_Rb + laser_freq)) ...
        / sqrt(zRx*zRy);
    
    Wx_Na = sqrt(-4*dipoleTrap0_unscaled_Na / (mass_Na * dipole_waist_x^2));
    Wx_Rb = sqrt(-4*dipoleTrap0_unscaled_Rb / (mass_Rb * dipole_waist_x^2));
    Wy = sqrt(-4*dipoleTrap0_unscaled_Na / (mass_Na * dipole_waist_y^2));
    
    %% making Udp(0) dimensionless
    dipoleTrap0_Na = dipoleTrap0_unscaled_Na / (getphysconst('hbar') * trap_freq);
    dipoleTrap0_Rb = dipoleTrap0_unscaled_Rb / (getphysconst('hbar') * trap_freq);
    
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