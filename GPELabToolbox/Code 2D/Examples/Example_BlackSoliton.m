%%% This file is an example of how to use GPELab (FFT version)

%% GROUND STATE COMPUTATION WITH A ROTATING TERM


%-----------------------------------------------------------
% Setting the data
%-----------------------------------------------------------

%% Setting the method and geometry
Computation = 'Ground';
Ncomponents = 1;
Type = 'BESP';
Deltat = 1e-1;
Stop_time = [];
Stop_crit = {'Energy',1e-9};
Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
xmin = -20;
xmax = 20;
ymin = -20;
ymax = 20;
Nx = 2^8+1;
Ny = 2^8+1;
Geometry2D = Geometry2D_Var2d(xmin,xmax,ymin,ymax,Nx,Ny);

%% Setting the physical problem
Delta = 0.5;
Beta = 10000;
Physics2D = Physics2D_Var2d(Method, Delta, Beta);
Physics2D = Dispersion_Var2d(Method, Physics2D);
Physics2D = Nonlinearity_Var2d(Method, Physics2D);

%% Setting the initial data
InitialData_Choice = 1;
Phi_0 = InitialData_Var2d(Method, Geometry2D, Physics2D, InitialData_Choice);
Phi_0{1} = Phi_0{1}.*(Geometry2D.X - 1i*Geometry2D.Y);


%% Setting informations and outputs
Outputs = OutputsINI_Var2d(Method);
Printing = 1;
Evo = 15;
Draw = 1;
Print = Print_Var2d(Printing,Evo,Draw);

%-----------------------------------------------------------
% Launching computation
%-----------------------------------------------------------

[Phi_1, Outputs] = GPELab2d(Phi_0,Method,Geometry2D,Physics2D,Outputs,[],Print);

%% Setting the method and geometry
Computation = 'Dynamic';
Ncomponents = 1;
Type = 'Relaxation';
Deltat = 1e-3;
Stop_time = 2;
Stop_crit = [];
Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);

%% Setting the initial data
eta = pi/2;
d = 0.4;
X_0 = 4;
Phi_1{1} = exp(1i*eta*tanh((Geometry2D.X-X_0)/d)).*Phi_1{1};

%% Setting informations and outputs
Save_solution = 1;
Save_evo = 5;
Outputs = OutputsINI_Var2d(Method,Save_evo,Save_solution);
Printing = 1;
Evo = 15;
Draw = 1;
Print = Print_Var2d(Printing,Evo,Draw);

%-----------------------------------------------------------
% Launching simulation
%-----------------------------------------------------------

[Phi, Outputs] = GPELab2d(Phi_1,Method,Geometry2D,Physics2D,Outputs,[],Print);