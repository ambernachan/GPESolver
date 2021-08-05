function [trap] = dipole_plus_quadratictrap(parameters, gx,gy,gz, X,Y,Z)

    dipole = dipoleTrap(parameters, 0,Y,Z);
    quadratic = quadraticTrap(parameters, gx,gy,gz, X,Y,Z);
    
    trap = cell(3,3);
    for n = 1:parameters.nComponents
        for m = 1:parameters.nComponents
            if n == m % diagonal elements
                trap{n,m} = quadratic{n,m} + dipole{n,m};
            else
                trap{n,m} = 0;
            end
        end
    end
    
end