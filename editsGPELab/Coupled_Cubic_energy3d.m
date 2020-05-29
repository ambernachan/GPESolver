function [CoupledCubicEnergy] = Coupled_Cubic_energy3d(Beta_coupled)
CoupledCubicEnergy = cell(2);
CoupledCubicEnergy{1,1} = @(Phi,X,Y,Z) 0.5 * Beta_coupled(1,1) * abs(Phi{1}).^2 + 0.5 * Beta_coupled(1,2) * abs(Phi{2}).^2;
CoupledCubicEnergy{2,2} = @(Phi,X,Y,Z) 0.5 * Beta_coupled(2,2) * abs(Phi{2}).^2 + 0.5 * Beta_coupled(2,1) * abs(Phi{1}).^2;
CoupledCubicEnergy{1,2} = @(Phi,X,Y,Z) 0;
CoupledCubicEnergy{2,1} = @(Phi,X,Y,Z) 0;
end