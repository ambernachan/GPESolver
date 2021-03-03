%{
...
%}

clear;
pause('on'); 

%% Determine interaction strength compared to kinetic energy

S = 100;

%% Saving info files

dimensions = 1;
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

Method = Method_Var1d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit, Max_iter);

% Saving workspace with relevant data for fitting
Method_ground = Method;
save(info.get_workspace_path('fittingdata'), 'Method_ground');
clear Method_ground;

%% Geometry1D

xmin = -15;
xmax = 15;
Nx = 2^9 + 1;
Geometry1D = Geometry1D_Var1d(xmin, xmax, Nx);

%% Physics1D

Delta = 0.5;
Beta = 4*pi*S;

Physics1D = Physics1D_Var1d(Method, Delta, Beta);
Physics1D = Dispersion_Var1d(Method, Physics1D);
Physics1D = Potential_Var1d(Method, Physics1D); % std quadratic potential
Physics1D = Nonlinearity_Var1d(Method, Physics1D); % std cubic nonlinearity


%% Defining a starting function Phi_0

InitialData_choice = 2; % TF initial data
w = (1+Beta/(2*pi))^(1/4); % interaction strength, w<1 for interactions; w=1 no interactions
s = 4; % size of the condensate for the trial function; s=1 is the expected result
X0 = 2;
gamma_x = 1 / (s*w)^2;

Phi_0 = InitialData_Var1d(Method, Geometry1D, Physics1D, InitialData_choice, X0, gamma_x);

%% Determining outputs
Outputs = OutputsINI_Var1d(Method);

%% Printing preliminary outputs
Printing = 1;
Evo = 15;
Draw = 0;
Print = Print_Var1d(Printing, Evo, Draw);

%% RUN THE SIMULATION to find the ground state

info.add_simulation_info(Geometry1D);
[Phi_1, Outputs] = GPELab1d(Phi_0, Method, Geometry1D, Physics1D, Outputs, [], Print);

%% Save the workspace & simulation info

info.add_result_info(Method, Outputs); 
save(info.get_workspace_path('groundstate'));

%% Draw & save solution

close all;
pause(2) % pauses the program for 2 seconds

Draw_solution1d(Phi_0, Method, Geometry1D, Figure_Var1d());

info.save_figure(1, 'initialdata', 'psi_sq');
info.save_figure(2, 'initialdata', 'angle');

Draw_solution1d(Phi_1, Method, Geometry1D, Figure_Var1d());

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

Method = Method_Var1d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);

% Saving workspace with relevant data for fitting
Method_dynamical = Method;
save(info.get_workspace_path('fittingdata'), 'Method_dynamical');
clear Method_dynamical;

%% We keep Geometry1D as-is

%% We keep Physics1D as-is

%% Determining outputs

Save_solution = 1;
Outputs = OutputsINI_Var1d(Method, Save_solution);

%% Printing preliminary outputs
Printing = 1;
Evo = 10;
Draw = 0;
Print = Print_Var1d(Printing, Evo, Draw);

%% RUN THE SIMULATION to find the ground state

[Phi, Outputs] = GPELab1d(Phi_1, Method, Geometry1D, Physics1D, Outputs, [], Print);

%% Save the workspace & simulation info

info.add_result_info(Method, Outputs);
info.finish_info();
save(info.get_workspace_path('dynamics'));

%% Draw & save solution

close all;
pause(10) % pauses the program for 10 seconds, program errors in next line, 
% but not when i put a breakpoint just before

Draw_solution1d(Phi, Method, Geometry1D, Figure_Var1d());

info.save_figure(1, 'dynamics', 'psi_sq');
info.save_figure(2, 'dynamics', 'angle');

%% Save PhiData structures for fitting w/o the whole workspace

phi_dyn = PhiData(Phi, Geometry1D);
phi_ground = PhiData(Phi_1, Geometry1D);
phi_input = PhiData(Phi_0, Geometry1D);

% Saving workspace with relevant data for fitting
save(info.get_workspace_path('fittingdata'), ... 
    'phi_dyn', 'phi_ground', 'phi_input', 'info', ... % necessary data
    'S', 'w', 'Beta', 'Method', ... % additional data
    '-append'); % to not overwrite Method_ground

%% end
