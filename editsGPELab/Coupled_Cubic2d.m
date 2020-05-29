function [CoupledCubicNonlinearity] = Coupled_Cubic2d(Method, Beta_coupled)

CoupledCubicNonlinearity = cell(2);
CoupledCubicNonlinearity{1,1} = @(Phi,X,Y) Beta_coupled(1,1)*abs(Phi{1}).^2 + Beta_coupled(1,2)*abs(Phi{2}).^2;
CoupledCubicNonlinearity{2,2} = @(Phi,X,Y) Beta_coupled(2,2)*abs(Phi{2}).^2 + Beta_coupled(2,1)*abs(Phi{1}).^2;
CoupledCubicNonlinearity{1,2} = @(Phi,X,Y) 0;
CoupledCubicNonlinearity{2,1} = @(Phi,X,Y) 0;
end