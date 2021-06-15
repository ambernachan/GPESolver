function [rho] = Population(Method, Geometry3D, Phi, X, Y, Z, FFTX, FFTY, FFTZ, m)
    
    Global_L2norm = 0;
    for n = 1:Method.Ncomponents
        Global_L2norm = Global_L2norm + L2_norm3d(Phi{n},Geometry3D)^2; % Computing the norm of each wave function
    end
    for n = 1:Method.Ncomponents
        Phi{n} = Phi{n} / sqrt(Global_L2norm); % Computing the norm of each wave function
    end
    
    TotalPhi = sum(sum(sum(abs(Phi{1})))).^2+sum(sum(sum(abs(Phi{2})))).^2+sum(sum(sum(abs(Phi{3})))).^2;
    rho_up = sum(sum(sum(abs(Phi{1})))).^2 / TotalPhi;
    rho_zero = sum(sum(sum(abs(Phi{2})))).^2 / TotalPhi;
    rho_down = sum(sum(sum(abs(Phi{3})))).^2 / TotalPhi;
    
    if m == 1
        rho = rho_up;
    elseif m == 0
        rho = rho_zero;
    elseif m == -1
        rho = rho_down;
    else
        error('Please specify the hyperfine projection mF=1,0,-1 as a final argument.')
    end
    
end