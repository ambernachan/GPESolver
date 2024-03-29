function [CoupledSpin1NonlinEnergy] = Coupled_CubicEnergy1d_spin1(Betan, Betas, params)
    
%     if ~(nargin > 2)
%         p = @(z) 0;
%         q = @(z) 0;
%     else
%         [p,q] = getVariableMagneticFieldPars(params.Bz, params.Bmin, params.trapfreq, params.Ehfs, params.boxlimits(1));
%     end
    
    % In the spin-1 manifold, we have 3 components of the wavefunction.
    % Phi{1} represents mF=1, Phi{2} reps mF=0 and Phi{3} reps mF=-1.
    CoupledSpin1NonlinEnergy = zerocell(3); % creates 3x3 cell of zeros
    
    % diagonal terms
    CoupledSpin1NonlinEnergy{1,1} = @(Phi,X) 0.5 * Betan * ...
        ( abs(Phi{1}).^2 + abs(Phi{2}).^2 + abs(Phi{3}).^2 ) + ...
        0.5 * Betas * ( abs(Phi{1}).^2 + abs(Phi{2}).^2 - abs(Phi{3}).^2 );
%         + (-p(X)+q(X)) .* abs(Phi{1}).^2;
%         + (q-p) * sum(sum(sum(abs(Phi{1}).^2)));
    CoupledSpin1NonlinEnergy{2,2} = @(Phi,X) 0.5 * Betan * ...
        ( abs(Phi{1}).^2 + abs(Phi{2}).^2 + abs(Phi{3}).^2 ) + ...
        0.5 * Betas * ( abs(Phi{1}).^2 + abs(Phi{3}).^2 );
    CoupledSpin1NonlinEnergy{3,3} = @(Phi,X) 0.5 * Betan * ...
        ( abs(Phi{1}).^2 + abs(Phi{2}).^2 + abs(Phi{3}).^2 ) + ...
        0.5 * Betas * ( abs(Phi{3}).^2 + abs(Phi{2}).^2 - abs(Phi{1} ).^2);
%         + (p(X)+q(X)) .* abs(Phi{3}).^2;
%         + (p+q) * sum(sum(sum(abs(Phi{3}).^2)));
    
    % spin-mixing terms
    % note that B = A' is the complex conjugate transpose of A
% %     CoupledSpin1NonlinEnergy{1,2} = @(Phi,X) Betas * Phi{3}' .* Phi{2};
% %     CoupledSpin1NonlinEnergy{2,1} = @(Phi,X) Betas * Phi{2}' .* Phi{3};
% %     CoupledSpin1NonlinEnergy{2,3} = @(Phi,X) Betas * Phi{2}' .* Phi{1};
% %     CoupledSpin1NonlinEnergy{3,2} = @(Phi,X) Betas * Phi{1}' .* Phi{2};

%     CoupledSpin1NonlinEnergy{1,2} = @(Phi,X) 0.5 * Betas * real(conj(Phi{3}) .* Phi{2});
%     CoupledSpin1NonlinEnergy{2,1} = @(Phi,X) 0.5 * Betas * real(conj(Phi{2}) .* Phi{3});
%     CoupledSpin1NonlinEnergy{2,3} = @(Phi,X) 0.5 * Betas * real(conj(Phi{2}) .* Phi{1});
%     CoupledSpin1NonlinEnergy{3,2} = @(Phi,X) 0.5 * Betas * real(conj(Phi{1}) .* Phi{2});    
    CoupledSpin1NonlinEnergy{1,2} = @(Phi,X) 0.5 * Betas * (conj(Phi{3}) .* Phi{2});
    CoupledSpin1NonlinEnergy{2,1} = @(Phi,X) 0.5 * Betas * (conj(Phi{2}) .* Phi{3});
    CoupledSpin1NonlinEnergy{2,3} = @(Phi,X) 0.5 * Betas * (conj(Phi{2}) .* Phi{1});
    CoupledSpin1NonlinEnergy{3,2} = @(Phi,X) 0.5 * Betas * (conj(Phi{1}) .* Phi{2});
    
    CoupledSpin1NonlinEnergy{3,1} = @(Phi,X) 0;
    CoupledSpin1NonlinEnergy{1,3} = @(Phi,X) 0;
    
    
end