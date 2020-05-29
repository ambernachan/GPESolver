clear all;
%%
Computation = 'Ground';
%Ncomponents = 1; % example case in GPE paper #1
Ncomponents = 2;
Type = 'BESP';
Deltat = 5e-2; % example case in GPE paper #1
%Deltat = 1e-5;
Stop_time = [];
Stop_crit = {'MaxNorm', 1e-6}; % example case in GPE paper #1
%Stop_crit = {'MaxNorm', 1e-5};

Method = Method_Var3d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);

%%

xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;
zmin = -15;
zmax = 15;

Nx = 2^7 + 1;
Ny = 2^7 + 1;
Nz = 2^7 + 1;

Geometry3D = Geometry3D_Var3d(xmin,xmax, ymin,ymax, zmin,zmax, Nx,Ny,Nz);
%Geometry2D = Geometry2D_Var2d(xmin,xmax, ymin,ymax, Nx,Ny);

%%
Delta = 0.5;
%Beta = 500;
Beta = 1000;
Beta_coupled = [1,2; 2,1];
%Omega = 0.1;
Omega = [0,0,0.7];
%Kappa = 1.75;
gammax = 1;
gammay = 1;
gammaz = 0.5;

Physics3D = Physics3D_Var3d(Method,Delta,Beta,Omega);

Physics3D = Dispersion_Var3d(Method, Physics3D); % example case in GPE paper #1
%Physics3D = Dispersion_Var3d(Method,Physics3D,Dispersion_SpinOrbit3d(Kappa));

% example case in GPE paper #1 :
Physics3D = Potential_Var3d(Method, Physics3D, @(X,Y,Z) quadratic_potential3d(gammax,gammay,gammaz, X,Y,Z));
%Physics3D = Potential_Var3d(Method, Physics3D);

Physics3D = Gradientx_Var3d(Method,Physics3D);
Physics3D = Gradienty_Var3d(Method,Physics3D);
Physics3D = Gradientz_Var3d(Method,Physics3D);

Physics3D = Nonlinearity_Var3d(Method,Physics3D,Coupled_Cubic3d(Beta_coupled),[],Coupled_Cubic_energy3d(Beta_coupled));
%Physics3D = Nonlinearity_Var3d(Method,Physics3D); % example case in GPE paper #1

%%

InitialData_choice = 2; % Thomas-Fermi approximation
Phi_0 = InitialData_Var3d(Method,Geometry3D,Physics3D,InitialData_choice);

Outputs = OutputsINI_Var3d(Method);
Printing = 1;
Evo = 15;
Draw = 1;
Print = Print_Var3d(Printing,Evo,Draw);

[Phi,Outputs] = GPELab3d(Phi_0,Method,Geometry3D,Physics3D,Outputs,[],Print);

%% end