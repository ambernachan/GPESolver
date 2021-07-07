function [fname] = createTextFileName(info, geometry, method, iterations)

    if strcmp(info.params.atom, 'Na')
        str_atom = '23Na';
    elseif strcmp(info.params.atom, 'Rb')
        str_atom = '87Rb';
    else
        str_atom = info.params.atom;
    end
    
    str_dx = sprintf('%.3g', geometry.dx);
    str_dt = sprintf('%.3g', method.Deltat);
    
    str_its = sprintf('%d', iterations);
    
    fname = [str_atom '_dx=' str_dx  '_dt=' str_dt '_its=' str_its];
    
end
    