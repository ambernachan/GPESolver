function [M] = Directional_Magnetization(Method, Geometry3D, Phi, direction, X, Y, Z, FFTX, FFTY, FFTZ)
    
    F = findExpF(Method, Geometry3D, Phi);
    
    Phi = normalize_global(Method, Geometry3D, Phi);
    TotalPhi = sum(sum(sum(abs(Phi{1}.^2))))+sum(sum(sum(abs(Phi{2}.^2))))+sum(sum(sum(abs(Phi{3}.^2))));
    for n = 1:Method.Ncomponents
%         p{n} = sum(sum(sum(Phi{n})));
%         p{n} = sum(sum(sum(abs(Phi{n}))));
        p{n} = sqrt(sum(sum(sum(abs(Phi{n}).^2)))); % absolute value of phi component
%         p1{n} = sqrt(sum(sum(sum(real(Phi{n}).^2))));
%         p2{n} = sqrt(sum(sum(sum(imag(Phi{n}).^2))));
%         P{n} = p1{n} + 1i * p2{n};
        P{n} = Phi{n} ./ sqrt(TotalPhi);
        Ph{n} = angle(Phi{n}); % phase of phi component
    end
    phi = [p{1}; p{2}; p{3}] ./ sqrt(TotalPhi);
    PHI = [P{1}; P{2}; P{3}] ./ sqrt(TotalPhi);
    
    Fx = (1/sqrt(2)) * [0 1 0; 1 0 1; 0 1 0];
    Fy = (1i/sqrt(2)) * [0 -1 0; 1 0 -1; 0 1 0];
    Fz = [1 0 0; 0 0 0; 0 0 -1];
    
    Mx = phi' * Fx * phi; 
    My = phi' * Fy * phi;
    Mz = phi' * Fz * phi;
    
    M = [{Mx}, {My}, {Mz}];
    
    % these small imaginary numbers are generally calculation errors;
    % remove them for clarity
    for i = 1:length(M)
        if imag(M{i}) < 1e-10
            M{i} = real(M{i});
        end
    end
%     
%     if strcmp(direction, 'x')
%         M = M{1};
%     elseif strcmp(direction, 'y')
%         M = M{2};
%     elseif strcmp(direction, 'z')
%         M = M{3};
%     else
%         error('No direction specified (fourth argument)')
%     end
%     
% %     for n = 1:3
% % %         Phi{n} = Phi{n} ./ sqrt(TotalPhi);
% %         abscoef{n} = abs(Phi{n});
% %     end
% %     phasediff12 = Ph{1} - Ph{2}; % (theta+ - theta0)
% %     phasediff23 = Ph{2} - Ph{3}; % (theta0 - theta-)
% %     
% %     MX = sqrt(2) * abscoef{2} .* ( abscoef{1} .* cos(phasediff12) + abscoef{3} .* cos(phasediff23));
% %     MY = -sqrt(2) * abscoef{2} .* ( abscoef{1} .* sin(phasediff12) + abscoef{3} .* sin(phasediff23));
% %     MZ = abscoef{1}.^2 - abscoef{3}.^2;
% %     
% %     mX = sum(sum(sum(MX)));
% %     mY = sum(sum(sum(MY)));
% %     mZ = sum(sum(sum(MZ)));

    if strcmp(direction, 'Mx')
        M = M{1};
    elseif strcmp(direction, 'My')
        M = M{2};
    elseif strcmp(direction, 'Mz')
        M = M{3};
    elseif strcmp(direction, 'F2')
        M = (M{1} + M{2} + M{3})^2;
    elseif strcmp(direction, 'M2')
        M = (M{1}^2 + M{2}^2 + M{3}^2);
    elseif strcmp(direction, 'x')
        M = F{1};
    elseif strcmp(direction, 'y')
        M = F{2};
    elseif strcmp(direction, 'z')
        M = F{3};
    else
        error('No direction specified (fourth argument)')
    end
    
end