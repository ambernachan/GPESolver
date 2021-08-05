function [trap] = quadraticTrap(parameters, gx,gy,gz, X,Y,Z)
    
    trap = cell(3,3);
    for n = 1:parameters.nComponents
        for m = 1:parameters.nComponents
            if n == m % diagonal elements
                trap{n,m} = quadratic_potential3d(gx, gy, gz, X, Y, Z);
            else
                trap{n,m} = 0;
            end
        end
    end
    
end