%{
Dynamics of a BEC in 3d
%}

fprintf('Start %s\n', datestr(now,'HH:MM:SS.FFF'))
clear;

%%
%{
GROUND STATE SIMULATION TO GET STARTING FUNCTION PHI_0 FOR DYNAMIC SIMULATION
%}

%% Simulation methods

Computation = 'Dynamic';
Ncomponents = 1;
Type = 'Splitting';
Deltat = 1e-3;
Stop_time = 1;
Stop_crit = [];

Method = Method_Var3d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);

%% Geometry3D

xmin = -2;
xmax = 2;
ymin = -2;
ymax = 2;
zmin = -2;
zmax = 2;
Nx = 2^7 + 1;
Ny = 2^7 + 1;
Nz = 2^7 + 1;
Geometry3D = Geometry3D_Var3d(xmin, xmax, ymin, ymax, zmin, zmax, Nx, Ny, Nz);

%% Physics3D for BEC with default

Delta = 1; % I think this should be 0.5...
Beta = 1e-3;
Physics3D = Physics3D_Var3d(Method, Delta, Beta);
Physics3D = Potential_Var3d(Method, Physics3D);
Physics3D = Nonlinearity_Var3d(Method, Physics3D);

%% Defining a starting function Phi_0

InitialData_choice = 2; % initialization by Thomas-Fermi approximation
Phi_0 = InitialData_Var3d(Method, Geometry3D, Physics3D, InitialData_choice);

%% Determining outputs
Outputs_iterations = 10;
Outputs_save = 0;
Outputs = OutputsINI_Var3d(Method, Outputs_iterations, Outputs_save);

%% Printing preliminary outputs
Printing = 1;
Evo = 15;
Draw = 1;
Print = Print_Var3d(Printing, Evo, Draw);

%% RUN THE SIMULATION to find the ground state

[Phi_1, Outputs] = GPELab3d(Phi_0, Method, Geometry3D, Physics3D, Outputs, [], Print);

%% Save the workspace

filename = mfilename;
string = filename(1:21); % Returns 'GPEpaper#_Example#_#D'

save(strcat(string, '-workspace-groundstate'))
save(strcat(string, '-outputs-groundstate'), 'Phi_1', 'Outputs')

%%
%{
DYNAMICAL SIMULATION
%}

%% We keep Mehod as-is

%%% Simulation methods
%
%Computation = 'Dynamic';
%Ncomponents = 1;
%Type = 'Relaxation';
%Deltat = 1e-3;
%Stop_time = 5;
%Stop_crit = [];
%
%Method = Method_Var3d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);

%% We keep Geometry3D as-is

%% Physics3D: we add a random phase to the initial state

A = 1;
d = 0.5;
Random_Phase = Stationary_Gaussian_Field3d(Geometry3D, ...
    @(X,Y,Z) A*exp(-(X.^2+Y.^2+Z.^2)/(2*d^2)));
Phi_1{1} = exp(-2i*pi*Random_Phase);

%% Determining outputs, including a user-defined function computing
%% the mean momentum of the BEC

Save_solution = 0;
Outputs_iterations = 10;

% Computing the mean momentum < \vec(p)_j >_{j=1,2,3} of the BEC
Output_function{1} = @(Phi,X,Y,Z,FFTX,FFTY,FFTZ) Geometry3D.dx*Geometry3D.dy*...
    Geometry3D.dz*sum( sum( sum( ifftn(FFTX.*fftn(Phi)).*conj(Phi) ) ) );
Output_function{2} = @(Phi,X,Y,Z,FFTX,FFTY,FFTZ) Geometry3D.dx*Geometry3D.dy*...
    Geometry3D.dz*sum( sum( sum( ifftn(FFTY.*fftn(Phi)).*conj(Phi) ) ) );
Output_function{3} = @(Phi,X,Y,Z,FFTX,FFTY,FFTZ) Geometry3D.dx*Geometry3D.dy*...
    Geometry3D.dz*sum( sum( sum( ifftn(FFTZ.*fftn(Phi)).*conj(Phi) ) ) );
Output_name{1} = 'BEC Momentum X';
Output_name{2} = 'BEC Momentum Y';
Output_name{3} = 'BEC Momentum Z';

Outputs = OutputsINI_Var3d(Method, Outputs_iterations, Save_solution, ...
    Output_function, Output_name);

%% Printing preliminary outputs
Printing = 1;
Evo = 10;
Draw = 1;
Print = Print_Var3d(Printing, Evo, Draw);

%% RUN THE SIMULATION to find the ground state

[Phi, Outputs] = GPELab3d(Phi_1, Method, Geometry3D, Physics3D, Outputs, [], Print);
fprintf('End %s\n', datestr(now,'HH:MM:SS.FFF'))

%% Save the workspace
save(strcat(string, '-workspace-dynamics'))
save(strcat(string, '-outputs-dynamics'), 'Phi', 'Outputs')

%% Drawing 1e-6 isovalues of the modulus with no transparency (alpha=1)

View = 3;
Isovalue = 1e-6;
Aspect = 1;
Figure = Figure_Var3d(View, Isovalue, Aspect);

Draw_Solution3d(Phi, Method, Geometry3D, Figure);

%% end
