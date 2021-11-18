function [x, labelx] = makedynamictimeline(x, info)

    x = x / info.params.trapfreq; % time in seconds

    if max(x) < 10^(-3) % max time < ms
        x = x * 10^6;
        labelx = 'time (\mus)';
    elseif (10^(-3) < max(x)) && (max(x) < 1) % max time > ms but < s
        x = x * 10^3;
        labelx = 'time (ms)';
    elseif max(x) > 1 % max time > s
        labelx = 'time (s)';
    end
end