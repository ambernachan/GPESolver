%{
...
%}

clear;

%% Determine interaction strength compared to kinetic energy

S = 0.01;

%% Saving info files
currentdir = pwd;
%workingFile = mfilename;
%fnameBase = workingFile(1:16); % = 'trial_Gaussian2D'
fnameBase = 'trial_Gaussian2D';

fnameInfo = strcat(currentdir, '/xOutputs/', fnameBase, '-info.txt'); % infofile name
fInfo = fopen(fnameInfo, 'w'); % create new file with simulation info
tStart = tic; % starting time for timer
CPUtimeAtStart = cputime; % CPU start time

fprintf( fInfo, 'Start: %s\n', datestr(now, 'dd mmm yy @ HH:MM:SS')); % print start time

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

%% Geometry2D

xmin = -5;
xmax = 5;
ymin = -5;
ymax = 5;
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

InitialData_choice = 1; % Gaussian initial data
w = (1 + Beta / (2*pi) )^(1/4); % interaction strength, w<1 for interactions; w=1 no interactions
s = 4; % size of the condensate for the trial function; s=1 is the expected result
X0 = 2;
Y0 = 2;
gamma_x = 1 / (s*w^2);
gamma_y = 1 / (s*w^2);

Phi_0 = InitialData_Var2d(Method, Geometry2D, Physics2D, InitialData_choice, X0, Y0, gamma_x, gamma_y);

%% Printing interaction strength

fprintf( fInfo, 'S =\t %f', S); % print interaction strength S
fprintf( fInfo, 'Beta =\t %f', Beta); % print interaction parameter Beta
fprintf( fInfo, 'w =\t %f', w); % print Gaussian parameter w

%% Determining outputs
Outputs = OutputsINI_Var2d(Method);

%% Printing preliminary outputs
Printing = 1;
Evo = 15;
Draw = 0;
Print = Print_Var2d(Printing, Evo, Draw);

%% RUN THE SIMULATION to find the ground state

[Phi_1, Outputs] = GPELab2d(Phi_0, Method, Geometry2D, Physics2D, Outputs, [], Print);

%% Save the workspace & simulation info
endTime = now;
tElapsed = toc(tStart); % save time elapsed
dimensions = 2; % get dimensionality from file name

% save information about final simulation iteration in info file
add_info( fInfo, Outputs, Method, dimensions, tElapsed, CPUtimeAtStart, cputime); 
fprintf( fInfo, 'Nx:\t%d\n', Nx);
fprintf( fInfo, 'Ny:\t%d\n', Ny);
fprintf( fInfo, 'End ground state sims: %s\n', datestr(endTime, 'dd mmm yy @ HH:MM:SS')); % end time

% save workspace to workspace folder
fnameWspace = strcat('xWorkspace-data/', fnameBase); % = 'GPEpaper#_Example#_#D'
save(strcat(fnameWspace, '-workspace-groundstate') )
save(strcat(fnameWspace, '-outputs-groundstate') , 'Phi_1', 'Outputs')

% (print) starting time for dynamical simulation
tStartDyn = tic; % starting time for dynamical sims timer
CPUtimeAtStartDyn = cputime; % CPU start time
fprintf( fInfo, 'Start dynamics: %s\n', datestr(now, 'dd mmm yy @ HH:MM:SS'));

%% Draw & save solution

fnameFigspace = strcat('xFigures/', fnameBase);

Draw_solution2d(Phi_0, Method, Geometry2D, Figure_Var2d());
savefig(1, strcat(fnameFigspace, '-psi_sq-testfunction'))
savefig(2, strcat(fnameFigspace, '-angle-testfunction'))

Draw_solution2d(Phi_1, Method, Geometry2D, Figure_Var2d());
savefig(1, strcat(fnameFigspace, '-psi_sq-groundstate'))
savefig(2, strcat(fnameFigspace, '-angle-groundstate'))

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
Draw = 1;
Print = Print_Var2d(Printing, Evo, Draw);

%% RUN THE SIMULATION to find the ground state

[Phi, Outputs] = GPELab2d(Phi_1, Method, Geometry2D, Physics2D, Outputs, [], Print);

%% Save the workspace & simulation info
endTime = now;
tElapsed = toc(tStart); % save total time elapsed
tElapsedDyn = toc(tStartDyn); % save time elapsed during dynamical sims

% save information about final simulation iteration in info file
add_info( fInfo, Outputs, Method, dimensions, tElapsedDyn, CPUtimeAtStartDyn, cputime); 

% extra time information on full simulation (ground+dynamics)
fprintf( fInfo, '-------------------------------------------\n');
fprintf( fInfo, 'Total CPU time:\t%8.2f\n', cputime - CPUtimeAtStart);
fprintf( fInfo, 'Total elapsed time:\t' );
fprintf( fInfo, print_time(tElapsed) );
fprintf( fInfo, '-------------------------------------------\n');

fprintf( fInfo, 'End: %s\n', datestr(endTime, 'dd mmm yy @ HH:MM:SS')); % end time
fclose( fInfo );

% save workspace to workspace folder
save(strcat(fnameWspace, '-workspace-dynamics') )
save(strcat(fnameWspace, '-outputs-dynamics') , 'Phi', 'Outputs')

%% Draw & save solution

Draw_solution2d(Phi, Method, Geometry2D, Figure_Var2d());
savefig(1, strcat(fnameFigspace, '-psi_sq-dynamics'))
savefig(2, strcat(fnameFigspace, '-angle-dynamics'))

%% expected solution

%expected = 1;
%X0 = 0;
%Y0 = 0;
%gamma_x = 1;
%gamma_y = 1;

%Phi_exp = InitialData_Var2d(Method, Geometry2D, Physics2D, expected, X0, Y0, gamma_x, gamma_y);

%Draw_solution2d(Phi_exp, Method, Geometry2D, Figure_Var2d());

%Draw_solution2d(cellfun(@minus, Phi_exp, Phi_1, 'UniformOutput',false), Method, Geometry2D, Figure_Var2d());

%% end
