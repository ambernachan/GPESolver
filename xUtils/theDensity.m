%% Computation of the root mean square of a two functions
%% INPUTS:
%%          Phi: functions whose density is computed (matrix)
%%          geometry: Structure containing variables concerning the geometry of the problem in 3D (structure) (see geometry_Var3d.m)
%% OUTPUT:
%%          density: density of phi (double)
%% FUNCTIONS USED:
%%          L2_norm3d: To compute the root mean square (line 11)

function [density] = theDensity(Phi, geometry)
    
    gridgap = geometry.dx;
    if isfield(geometry, 'dy')
        gridgap = gridgap * geometry.dy;
        if isfield(geometry, 'dz')
            gridgap = gridgap * geometry.dz;
        end
    end
    
    sqphis = 0;
    if iscell(Phi)
        for n = 1:length(Phi)
            sqphis = sqphis + abs(Phi{n}).^2;
        end
    else
        sqphis = abs(Phi).^2;
    end
    density = (gridgap)*sum(sqphis.^2, 'all'); % Computation of the L2 norm
    
end