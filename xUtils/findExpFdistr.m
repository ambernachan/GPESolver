% gives normalized Fx,Fy,Fz distribution in (1x3) cell of (Nx*Ny*Nz) arrays
function [expectedFdistr] = findExpFdistr(Method, Geometry3D, Phi)
    
    Phi = normalize_global(Method, Geometry3D, Phi);
    TotalPhi = sum(sum(sum(abs(Phi{1}.^2))))+sum(sum(sum(abs(Phi{2}.^2))))+sum(sum(sum(abs(Phi{3}.^2))));
    expectedFdistr = cell(1,3);
    
    Fx = sqrt(2) * abs(Phi{2}) .* ( ...
        abs(Phi{1}) .* cos( angle(Phi{1})-angle(Phi{2}) ) + ...
        abs(Phi{3}) .* cos( angle(Phi{3})-angle(Phi{2}) ) ...
        );
    
    Fy = -sqrt(2) * abs(Phi{2}) .* ( ...
        abs(Phi{1}) .* sin( angle(Phi{1})-angle(Phi{2}) ) - ...
        abs(Phi{3}) .* sin( angle(Phi{3})-angle(Phi{2}) ) ...
        );
    
    Fz = abs(Phi{1}).^2 - abs(Phi{3}).^2;
    
    expectedFdistr{1} = Fx ./ TotalPhi;
    expectedFdistr{2} = Fy ./ TotalPhi;
    expectedFdistr{3} = Fz ./ TotalPhi;
    
end