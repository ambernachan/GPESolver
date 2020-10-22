%{
...
%}

clear;
pause('on'); 

%% Determine interaction strength compared to kinetic energy

S = 100;

%% Saving info files

dimensions = 3;
info = Info(name_from_filename(mfilename), dimensions);

%%
%{
GROUND STATE SIMULATION TO GET STARTING FUNCTION PHI_0 FOR DYNAMIC SIMULATION
%}

%% Simulation methods

Computation = 'Ground';
Ncomponents = 1;
Type = 'BESP';
Deltat = 0.05;
Stop_time = [];
Stop_crit = {'MaxNorm', 1e-4};
Max_iter = 6e4;

Method = Method_Var3d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit, Max_iter);

% Saving workspace with relevant data for fitting
Method_ground = Method;
save(info.get_workspace_path('fittingdata'), 'Method_ground');
clear Method_ground;

%% Geometry3D

xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;
zmin = -10;
zmax = 10;
Nx = 2^6 + 1;
Ny = 2^6 + 1;
Nz = 2^6 + 1;
Geometry3D = Geometry3D_Var3d(xmin, xmax, ymin, ymax, zmin, zmax, Nx, Ny, Nz);

%% Physics3D

Delta = 0.5;
Beta = 4*pi*S;
Omega = [0 0 0];
Physics3D = Physics3D_Var3d(Method, Delta, Beta, Omega);
Physics3D = Potential_Var3d(Method, Physics3D); % std quadratic potential
Physics3D = Nonlinearity_Var3d(Method, Physics3D); % std cubic nonlinearity
%Physics3D = Gradientx_Var3d(Method, Physics3D, @(x,y,z) -1i*Omega*y);(!!!)
%Physics3D = Gradienty_Var3d(Method, Physics3D, @(x,y,z) 1i*Omega*x); (!!!)

%% Defining a starting function Phi_0

InitialData_choice = 2; % TF initial data
w = (1+Beta/(2*pi))^(1/4); % interaction strength, w<1 for interactions; w=1 no interactions
s = 4; % size of the condensate for the trial function; s=1 is the expected result
X0 = 2;
Y0 = 2;
Z0 = 2;
%gamma_x = 1 / (s*w)^2;
%gamma_y = 1 / (s*w)^2;
%gamma_z = 1 / (s*w)^2;
gamma_x = 1;
gamma_y = 1;
gamma_z = 1;

Phi_0 = InitialData_Var3d(Method, Geometry3D, Physics3D, InitialData_choice, X0, Y0, Z0, gamma_x, gamma_y, gamma_z);

%% Determining outputs
Outputs = OutputsINI_Var3d(Method);
Evo = 15; Save = 1;
%Outputs = OutputsINI_Var3d(Method, Evo, Save, @(Phi) output_phisquared(Phi), 'phi_squared');
Outputs = OutputsINI_Var3d(Method, Evo, Save);

%% Printing preliminary outputs
Printing = 1;
Evo = 15;
Draw = 0;
Print = Print_Var3d(Printing, Evo, Draw);

%% RUN THE SIMULATION to find the ground state

info.add_simulation_info(Geometry3D);
[Phi_1, Outputs] = GPELab3d(Phi_0, Method, Geometry3D, Physics3D, Outputs, [], Print);

MakeVideo3d(Method,Geometry3D,Outputs); % std function |phi|^2

%% Save the workspace & simulation info

info.add_result_info(Method, Outputs); 
save(info.get_workspace_path('groundstate'));

%% Draw & save solution

close all;
pause(2) % pauses the program for 2 seconds

Draw_solution3d(Phi_0, Method, Geometry3D, Figure_Var3d());

info.save_figure(1, 'initialdata', 'psi_sq');
info.save_figure(2, 'initialdata', 'angle');

Draw_solution3d(Phi_1, Method, Geometry3D, Figure_Var3d());

info.save_figure(1, 'groundstate', 'psi_sq');
info.save_figure(2, 'groundstate', 'angle');

%%
%{
DYNAMICAL SIMULATION
%}

%% Simulation methods

Computation = 'Dynamic';
Ncomponents = 1;
Type = 'Splitting';
Deltat = 1e-4;
Stop_time = 1; % = 1e4 iterations max
Stop_time = [];
%Stop_crit = {'MaxNorm', 1e-4};
Stop_crit = {'Energy', 1e-8};

Method = Method_Var3d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);

% Saving workspace with relevant data for fitting
Method_dynamical = Method;
save(info.get_workspace_path('fittingdata'), 'Method_dynamical');
clear Method_dynamical;

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

%% Save the workspace & simulation info

info.add_result_info(Method, Outputs);
info.finish_info();
save(info.get_workspace_path('dynamics'));

%% Draw & save solution

close all;
pause(10) % pauses the program for 10 seconds, program errors in next line, 
% but not when i put a breakpoint just before

Draw_solution3d(Phi, Method, Geometry3D, Figure_Var3d());

info.save_figure(1, 'dynamics', 'psi_sq');
info.save_figure(2, 'dynamics', 'angle');

%% Save PhiData structures for fitting w/o the whole workspace

phi_dyn = PhiData(Phi, Geometry3D);
phi_ground = PhiData(Phi_1, Geometry3D);
phi_input = PhiData(Phi_0, Geometry3D);

% Saving workspace with relevant data for fitting
save(info.get_workspace_path('fittingdata'), ... 
    'phi_dyn', 'phi_ground', 'phi_input', 'info', ... % necessary data
    'S', 'w', 'Beta', 'Method', ... % additional data
    '-append'); % to not overwrite Method_ground

%% end
