%% Computation of global norm, returns Phi normalized over its components
%% INPUTS:
%%          Phi: Wave function whose L2 norm is computed (cell matrix)
%%          Geometry: As in struct from GPELab
%% OUTPUT:
%%          Phi: L2 norm of the function phi (double)
%% DEPENDENCIES:
%%          L2norm1d, L2norm2d, L2norm3d, respectively, for 1,2,3D
%% Normalization as:
%%          norm = Sum[ L2norm(Phi{n})^2 ], summed over n components
%%          Phi = Phi / sqrt(norm)

function [Phi] = normalize_global(Method, Geometry, Phi)
    
    nComponents = numel(Phi);
    
    Global_L2norm = 0;
    dims = 1;
    if isfield(Geometry, 'dy')
        dims = dims + 1;
        if isfield(Geometry, 'dz')
            dims = dims + 1;
        end
    end
%     norm_f_name = ['L2_norm' num2str(dims) 'd'];
    
    for n = 1:nComponents
        if dims == 1
            Global_L2norm = Global_L2norm + L2_norm1d(Phi{n},Geometry)^2; % Computing the norm of each wave function
        elseif dims == 2
            Global_L2norm = Global_L2norm + L2_norm2d(Phi{n},Geometry)^2; % Computing the norm of each wave function
        elseif dims == 3
            Global_L2norm = Global_L2norm + L2_norm3d(Phi{n},Geometry)^2; % Computing the norm of each wave function
        end
    end
    for n = 1:nComponents
        Phi{n} = Phi{n} / sqrt(Global_L2norm); % Computing the norm of each wave function
    end
    
end