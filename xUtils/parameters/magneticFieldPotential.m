function [potential] = magneticFieldPotential(parameters, Bz, Wmin, l)
    
%     potential = zeros(parameters.nComponents, parameters.nComponents);
%     potential = num2cell(potential);
%     
    [p, q] = getMagneticFieldPars(Bz, Wmin, parameters.Ehfs);
    % For component 1, mF=1 and Vmagn = (-p + q) in dimensionless notation
    % For component 3, mF=-1 and Vmagn = (p + q) in dimensionless notation
%     potential{1,1} = -p + q;
%     potential{3,3} = p + q;
    potential = - l * p + l^2 * q;
    
end