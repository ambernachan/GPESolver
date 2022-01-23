function [CoupledSpin1Nonlin] = Coupled_Cubic1d_spin1(Betan, Betas, params)
    
%     if ~(nargin > 2)
%         p = @(z) 0;
%         q = @(z) 0;
%     else
%         [p,q] = getVariableMagneticFieldPars(params.Bz, params.Bmin, params.trapfreq, params.Ehfs, params.boxlimits(1));
%     end
    
    % In the spin-1 manifold, we have 3 components of the wavefunction.
    % Phi{1} represents mF=1, Phi{2} reps mF=0 and Phi{3} reps mF=-1.
    CoupledSpin1Nonlin = zerocell(3); % creates 3x3 cell of @(Phi,X) 0's
    
    % diagonal terms
    CoupledSpin1Nonlin{1,1} = @(Phi,X) Betan * ...
        ( abs(Phi{1}).^2 + abs(Phi{2}).^2 + abs(Phi{3}).^2 ) + ...
        Betas * ( abs(Phi{1}).^2 + abs(Phi{2}).^2 - abs(Phi{3}).^2 );
%         - p(X) + q(X);
    CoupledSpin1Nonlin{2,2} = @(Phi,X) Betan * ...
        ( abs(Phi{1}).^2 + abs(Phi{2}).^2 + abs(Phi{3}).^2 ) + ...
        Betas * ( abs(Phi{1}).^2 + abs(Phi{3}).^2 );
    CoupledSpin1Nonlin{3,3} = @(Phi,X) Betan * ...
        ( abs(Phi{1}).^2 + abs(Phi{2}).^2 + abs(Phi{3}).^2 ) + ...
        Betas * ( abs(Phi{3}).^2 + abs(Phi{2}).^2 - abs(Phi{1} ).^2);
%         + p(X) + q(X);
    
    % spin-mixing terms
    % note that B = A' is the complex conjugate transpose of A
%     CoupledSpin1Nonlin{1,2} = @(Phi,X) Betas * Phi{3}' .* Phi{2};
%     CoupledSpin1Nonlin{2,1} = @(Phi,X) Betas * Phi{2}' .* Phi{3};
%     CoupledSpin1Nonlin{2,3} = @(Phi,X) Betas * Phi{2}' .* Phi{1};
%     CoupledSpin1Nonlin{3,2} = @(Phi,X) Betas * Phi{1}' .* Phi{2};
    CoupledSpin1Nonlin{1,2} = @(Phi,X) Betas * conj(Phi{3}) .* Phi{2};
    CoupledSpin1Nonlin{2,1} = @(Phi,X) Betas * conj(Phi{2}) .* Phi{3};
    CoupledSpin1Nonlin{2,3} = @(Phi,X) Betas * conj(Phi{2}) .* Phi{1};
    CoupledSpin1Nonlin{3,2} = @(Phi,X) Betas * conj(Phi{1}) .* Phi{2};
    
    CoupledSpin1Nonlin{3,1} = @(Phi,X) 0;
    CoupledSpin1Nonlin{1,3} = @(Phi,X) 0;
    
end