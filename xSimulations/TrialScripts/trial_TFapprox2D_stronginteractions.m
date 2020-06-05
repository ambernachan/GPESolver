%{
...
%}

clear;

%% Determine interaction strength compared to kinetic energy

S = 100;

%% Saving info files

dimensions = 2;
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
Stop_crit = {'MaxNorm', 1e-8};
Max_iter = 6e4;

Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit, Max_iter);

% Saving workspace with relevant data for fitting
Method_ground = Method;
save(info.get_workspace_path('fittingdata'), 'Method_ground');
clear Method_ground;

%% Geometry2D

xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;
Nx = 2^8 + 1;
Ny = 2^8 + 1;
Geometry2D = Geometry2D_Var2d(xmin, xmax, ymin, ymax, Nx, Ny);

%% Physics2D

Delta = 0.5;
Beta = 4*pi*S;
Omega = 0;
Physics2D = Physics2D_Var2d(Method, Delta, Beta, Omega);
Physics2D = Potential_Var2d(Method, Physics2D); % std quadratic potential
Physics2D = Nonlinearity_Var2d(Method, Physics2D); % std cubic nonlinearity
%Physics2D = Gradientx_Var2d(Method, Physics2D, @(x,y) -1i*Omega*y);
%Physics2D = Gradienty_Var2d(Method, Physics2D, @(x,y) 1i*Omega*x);

%% Defining a starting function Phi_0

InitialData_choice = 2; % Gaussian initial data
w = (1+Beta/(2*pi))^(1/4); % interaction strength, w<1 for interactions; w=1 no interactions
s = 4; % size of the condensate for the trial function; s=1 is the expected result
X0 = 2;
Y0 = 2;
gamma_x = 1 / (s*w)^2;
gamma_y = 1 / (s*w)^2;

Phi_0 = InitialData_Var2d(Method, Geometry2D, Physics2D, InitialData_choice, X0, Y0, gamma_x, gamma_y);

%% Determining outputs
Outputs = OutputsINI_Var2d(Method);

%% Printing preliminary outputs
Printing = 1;
Evo = 15;
Draw = 0;
Print = Print_Var2d(Printing, Evo, Draw);

%% RUN THE SIMULATION to find the ground state

info.add_simulation_info(Geometry2D);
[Phi_1, Outputs] = GPELab2d(Phi_0, Method, Geometry2D, Physics2D, Outputs, [], Print);

%% Save the workspace & simulation info

info.add_result_info(Method, Outputs); 
save(info.get_workspace_path('groundstate'));

%% Draw & save solution

Draw_solution2d(Phi_0, Method, Geometry2D, Figure_Var2d());

info.save_figure(1, 'initialdata', 'psi_sq');
info.save_figure(2, 'initialdata', 'angle');

Draw_solution2d(Phi_1, Method, Geometry2D, Figure_Var2d());

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
Stop_time = 1;
Stop_crit = [];

Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);

%% We keep Geometry2D as-is

%% We keep Physics2D as-is

%% Determining outputs

Save_solution = 1;
Outputs = OutputsINI_Var2d(Method, Save_solution);

%% Printing preliminary outputs
Printing = 1;
Evo = 10;
Draw = 0;
Print = Print_Var2d(Printing, Evo, Draw);

%% RUN THE SIMULATION to find the ground state

[Phi, Outputs] = GPELab2d(Phi_1, Method, Geometry2D, Physics2D, Outputs, [], Print);

%% Save the workspace & simulation info

info.add_result_info(Method, Outputs);
info.finish_info();
save(info.get_workspace_path('dynamics'));

%% Draw & save solution

Draw_solution2d(Phi, Method, Geometry2D, Figure_Var2d());

info.save_figure(1, 'dynamics', 'psi_sq');
info.save_figure(2, 'dynamics', 'angle');

%% Save PhiData structures for fitting w/o the whole workspace

phi_dyn = PhiData(Phi, Geometry2D);
phi_ground = PhiData(Phi_1, Geometry2D);
phi_input = PhiData(Phi_0, Geometry2D);

% Saving workspace with relevant data for fitting
save(info.get_workspace_path('fittingdata'), ... 
    'phi_dyn', 'phi_ground', 'phi_input', 'info', ... % necessary data
    'S', 'w', 'Beta', 'Method', ... % additional data
    '-append'); % to not overwrite Method_ground

%% end
