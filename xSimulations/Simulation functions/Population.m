function [rho] = Population(Method, geometry, Phi, X, Y, Z, FFTX, FFTY, FFTZ, m)
    
    Phi = normalize_global(Method, geometry, Phi);
    
    TotalPhi = sum(abs(Phi{1}.^2),'all')+sum(abs(Phi{2}.^2),'all')+sum(abs(Phi{3}.^2),'all');
    rho_up = sum(abs(Phi{1}.^2),'all') / TotalPhi;
    rho_zero = sum(abs(Phi{2}.^2),'all') / TotalPhi;
    rho_down = sum(abs(Phi{3}.^2),'all') / TotalPhi;
    
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