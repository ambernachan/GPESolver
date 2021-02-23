%{
Low-interactions system without rotations in 3D
%}

clear all;

%% Determine interaction strength compared to kinetic energy

S = 0; % parameter S is chi = N * a_s / a_ho, so Beta = 4*pi*chi


%%
%{
GROUND STATE SIMULATION TO GET STARTING FUNCTION PHI_0 FOR DYNAMIC SIMULATION
%}

%% Simulation methods

Computation = 'Ground';
Ncomponents = 1;
Type = 'BESP';
Deltat = 0.5;
Stop_time = [];
Stop_crit = {'MaxNorm', 1e-4};
Max_iter = 1e3;

Method = Method_Var3d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit, Max_iter);


%% Geometry3D

xmin = -5;
xmax = 5;
ymin = -5;
ymax = 5;
zmin = -5;
zmax = 5;
Nx = 2^7 + 1;
Ny = 2^7 + 1;
Nz = 2^7 + 1;

Geometry3D = Geometry3D_Var3d(xmin, xmax, ymin, ymax, zmin, zmax, Nx, Ny, Nz);

%% Physics3D

% Harmonic trap parameters
gx = 1;
gy = 1;
gz = 1;

Delta = 0.5;
Beta = 4*pi*S; % S is 'chi' in literature -> N * a_s / a_ho
Omega = 0;
Physics3D = Physics3D_Var3d(Method, Delta, Beta, Omega);
%Physics3D = Potential_Var3d(Method, Physics3D); % std quadratic potential
Physics3D = Potential_Var3d(Method, Physics3D, @(X,Y,Z) quadratic_potential3d(gx, gy, gz, X, Y, Z));
Physics3D = Nonlinearity_Var3d(Method, Physics3D); % std cubic nonlinearity

% no rotations
%Physics3D = Gradientx_Var3d(Method, Physics3D, @(x,y) -1i*Omega*y);
%Physics3D = Gradienty_Var3d(Method, Physics3D, @(x,y) 1i*Omega*x);

%% Defining a starting function Phi_0

InitialData_choice = 1; % Gaussian initial data

X0 = 0;
Y0 = 0;
Z0 = 0;
gamma_x = 1;
gamma_y = 1;
gamma_z = 1;

Phi_0 = InitialData_Var3d(Method, Geometry3D, Physics3D, InitialData_choice, X0, Y0, Z0, gamma_x, gamma_y, gamma_z);

%% Determining outputs

Outputs = OutputsINI_Var3d(Method);

%% Printing preliminary outputs

Printing = 1;
Evo = 15;
Draw = 0;
Print = Print_Var3d(Printing, Evo, Draw);

%% RUN THE SIMULATION to find the ground state

[Phi_1, Outputs] = GPELab3d(Phi_0, Method, Geometry3D, Physics3D, Outputs, [], Print);


%% Draw solution

close all;
pause(2) % pauses the program for 2 seconds

Draw_solution3d(Phi_1, Method, Geometry3D, Figure_Var3d());


%%
%{
DYNAMICAL SIMULATION
%}

%% Simulation methods

Computation = 'Dynamic';
Ncomponents = 1;
Type = 'Splitting';
Deltat = 1e-4;
Max_iter = 200; % testing value
Stop_crit = {'MaxNorm', 1e-4};

%Method = Method_Var3d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit, Max_iter);
Method = Method_Var3d(Computation, Ncomponents, Type, Deltat, [], Stop_crit, Max_iter);

%% We keep Geometry3D as-is

%% We keep Physics3D as-is

%% Determining outputs

Save_solution = 1;
Outputs = OutputsINI_Var3d(Method, Save_solution);

%% Printing preliminary outputs
Printing = 1;
Evo = 10;
Draw = 0;
Print = Print_Var3d(Printing, Evo, Draw);

%% RUN THE SIMULATION to find the ground state

[Phi, Outputs] = GPELab3d(Phi_1, Method, Geometry3D, Physics3D, Outputs, [], Print);

%% Draw & save solution

close all;
pause(2) % pauses the program for 2 seconds

Draw_solution3d(Phi, Method, Geometry3D, Figure_Var3d());

%% Save workspace
chistr = sprintf('chi=%.2f', S);
gammastr = sprintf('gamma=[%.1f,%.1f,%.1f]', gx, gy, gz);
gridstr = sprintf('grid[-%d,%d]', xmax, xmax);
workspacename = ['Gaussian3D_weakinteractions_MWE_' chistr '_' gammastr '_' gridstr '.mat']
save(workspacename)

%% end
