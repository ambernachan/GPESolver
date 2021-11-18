function [x, labelx] = makexaxisinmeters(x, info)

%     aho = sqrt( getphysconst('hbar') / (info.params.atom_mass * info.params.trapfreq) ); % characteristic h.o. length
    x = x * info.params.aho; % position in meters

    if max(x) < 1.1*10^(-3) % max position < mm
        x = x * 10^6;
        labelx = ' (\mum)';
    elseif (10^(-3) < max(x)) && (max(x) < 1.1) % max position > mm but < m
        x = x * 10^3;
        labelx = ' (mm)';
    elseif max(x) > 1 % max position > m
        labelx = ' (m)';
    end
end