function [Phi] = normalize_global(Method, Geometry, Phi)

    Global_L2norm = 0;
    for n = 1:Method.Ncomponents
        Global_L2norm = Global_L2norm + L2_norm3d(Phi{n},Geometry)^2; % Computing the norm of each wave function
    end
    for n = 1:Method.Ncomponents
        Phi{n} = Phi{n} / sqrt(Global_L2norm); % Computing the norm of each wave function
    end
    
end