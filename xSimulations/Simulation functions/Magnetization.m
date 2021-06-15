function [M] = Magnetization(Method, Geometry3D, Phi, X, Y, Z, FFTX, FFTY, FFTZ)
    
    Global_L2norm = 0;
    for n = 1:Method.Ncomponents
        Global_L2norm = Global_L2norm + L2_norm3d(Phi{n},Geometry3D)^2; % Computing the norm of each wave function
    end
    for n = 1:Method.Ncomponents
        Phi{n} = Phi{n} / sqrt(Global_L2norm); % Computing the norm of each wave function
    end
    
    TotalPhi = sum(sum(sum(abs(Phi{1})))).^2+sum(sum(sum(abs(Phi{2})))).^2+sum(sum(sum(abs(Phi{3})))).^2;
    Mup = sum(sum(sum(abs(Phi{1})))).^2;
    Mdown = sum(sum(sum(abs(Phi{3})))).^2;
    
    M = (Mup - Mdown) / TotalPhi; 
    
end