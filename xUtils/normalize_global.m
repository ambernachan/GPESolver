function [Phi] = normalize_global(Method, Geometry, Phi)

    Global_L2norm = 0;
    dims = 1;
    if isfield(Geometry, 'dy')
        dims = dims + 1;
        if isfield(Geometry, 'dz')
            dims = dims + 1;
        end
    end
%     norm_f_name = ['L2_norm' num2str(dims) 'd'];
    
    for n = 1:Method.Ncomponents
        if dims == 1
            Global_L2norm = Global_L2norm + L2_norm1d(Phi{n},Geometry)^2; % Computing the norm of each wave function
        elseif dims == 2
            Global_L2norm = Global_L2norm + L2_norm2d(Phi{n},Geometry)^2; % Computing the norm of each wave function
        elseif dims == 3
            Global_L2norm = Global_L2norm + L2_norm3d(Phi{n},Geometry)^2; % Computing the norm of each wave function
        end
    end
    for n = 1:Method.Ncomponents
        Phi{n} = Phi{n} / sqrt(Global_L2norm); % Computing the norm of each wave function
    end
    
end