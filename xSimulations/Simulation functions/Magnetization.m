function [M] = Magnetization(Method, Geometry3D, Phi, X, Y, Z, FFTX, FFTY, FFTZ)
    
    Phi = normalize_global(Method, Geometry3D, Phi);
    
    TotalPhi = sum(sum(sum(abs(Phi{1}.^2))))+sum(sum(sum(abs(Phi{2}.^2))))+sum(sum(sum(abs(Phi{3}.^2))));
    Mup = sum(sum(sum(abs(Phi{1}.^2))));
    Mdown = sum(sum(sum(abs(Phi{3}.^2))));
    
    M = (Mup - Mdown) / TotalPhi; 
    
end