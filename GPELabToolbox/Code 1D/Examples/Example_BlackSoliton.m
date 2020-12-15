%%% This file is an example of how to use GPELab (FFT version)

%% Ground state of a Gross-Pitaevskii equation with quadratic potential and cubic nonlinearity in 1D


clear all;
%-----------------------------------------------------------
% Setting the data
%-----------------------------------------------------------

%% Setting the method and geometry
Computation = 'Ground';
Ncomponents = 1;
Type = 'BESP';
Deltat = 1e-2;
Stop_time = [];
Stop_crit = {'MaxNorm',1e-9};
Method = Method_Var1d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
xmin = -15;
xmax = 15;
Nx = 2^11+1;
Geometry1D = Geometry1D_Var1d(xmin,xmax,Nx);

%% Setting the physical problem
Delta = 1;
Beta = 200;
Physics1D = Physics1D_Var1d(Method,Delta,Beta); 
Physics1D = Dispersion_Var1d(Method, Physics1D);
% Physics1D = Potential_Var1d(Method, Physics1D);
Physics1D = Nonlinearity_Var1d(Method, Physics1D); 

%% Setting the initial data
InitialData_Choice = 1;
%Phi_0 = InitialData_Var1d(Method, Geometry1D, Physics1D, InitialData_Choice);
Phi_0{1} = ones(Nx+1,1);

%% Setting informations and outputs
Solution_save = 1;
Outputs_iterations = 10;
Output_function{1} = @(Phi,X,FFTX) Geometry1D.dx*sum(X.*abs(Phi).^2);
Output_name{1} = 'Position of the soliton';
Outputs = OutputsINI_Var1d(Method,Outputs_iterations,Solution_save,Output_function,Output_name);
Printing = 1;
Evo = 100;
Draw = 1;
Print = Print_Var1d(Printing,Evo,Draw);

%-----------------------------------------------------------
% Launching simulation
%-----------------------------------------------------------

[Phi_1, Outputs] = GPELab1d(Phi_0,Method,Geometry1D,Physics1D,Outputs,[],Print);

%% Setting the initial data
L_h = 0.3;
c = sqrt(Beta/2);
v = 0.2*c;
Xi = sqrt(1-(v/c)^2);
X_0= 2;
X = Geometry1D.X ; 
Phi_1{1} = (Xi*tanh(Xi*X/L_h)+1i*(v/c)).*Phi_1{1};

%% Setting the method and geometry
Computation = 'Dynamic';
Ncomponents = 1;
Type = 'Relaxation';
Deltat = 1e-3;
Stop_time = 10;
Stop_crit = [];
Method = Method_Var1d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
xmin = -15;
xmax = 15;
Nx = 2^11+1;
Geometry1D = Geometry1D_Var1d(xmin,xmax,Nx);

%% Setting the physical problem
Delta = 1;
Beta = 200;
Physics1D = Physics1D_Var1d(Method,Delta,Beta);
Physics1D = Dispersion_Var1d(Method, Physics1D);
% Physics1D = Potential_Var1d(Method, Physics1D);
Physics1D = Nonlinearity_Var1d(Method, Physics1D); 

%% Setting informations and outputs
Solution_save = 1;
Outputs_iterations = 10;
Output_function{1} = @(Phi,X,FFTX) Geometry1D.dx*sum(X.*abs(Phi).^2);
Output_name{1} = 'Position of the soliton';
Outputs = OutputsINI_Var1d(Method,Outputs_iterations,Solution_save,Output_function,Output_name);

[Phi, Outputs] = GPELab1d(Phi_1,Method,Geometry1D,Physics1D,Outputs,[],Print);

Draw_Timesolution1d(Outputs,Method,Geometry1D,Figure_Var1d)