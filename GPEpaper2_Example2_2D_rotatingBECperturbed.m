%{
Dynamics of a rotating Bose-Einstein condensate perturbed by
a random Gaussian potential in 2D
%}

clear;

%% Saving info files

workingFile = mfilename;
fnameBase = workingFile(1:21); % = 'GPEpaper#_Example#_#D'

fnameInfo = strcat('xOutputs/', fnameBase, '-info.txt'); % infofile name
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
Deltat = 0.5;
Stop_time = [];
Stop_crit = {'MaxNorm', 1e-8};
Max_iter = 6e4;

Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit, Max_iter);

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
Beta = 1000;
Omega = 0.52;
Physics2D = Physics2D_Var2d(Method, Delta, Beta, Omega);
Physics2D = Potential_Var2d(Method, Physics2D, @(x,y) (1/2)*(x.^2+y.^2));
Physics2D = Nonlinearity_Var2d(Method, Physics2D, @(phi,x,y) abs(phi).^2);
Physics2D = Gradientx_Var2d(Method, Physics2D, @(x,y) -1i*Omega*y);
Physics2D = Gradienty_Var2d(Method, Physics2D, @(x,y) 1i*Omega*x);

%% Defining a starting function Phi_0

InitialData_choice = 2;
Phi_0 = InitialData_Var2d(Method, Geometry2D, Physics2D, InitialData_choice);

%% Determining outputs
Outputs = OutputsINI_Var2d(Method);

%% Printing preliminary outputs
Printing = 1;
Evo = 15;
Draw = 1;
Print = Print_Var2d(Printing, Evo, Draw);

%% RUN THE SIMULATION to find the ground state

[Phi_1, Outputs] = GPELab2d(Phi_0, Method, Geometry2D, Physics2D, Outputs, [], Print);

%% Save the workspace & simulation info
endTime = now;
tElapsed = toc(tStart); % save time elapsed
dimensions = str2double(fnameBase(20)); % get dimensionality from file name

% save information about final simulation iteration in info file
add_info( fInfo, Outputs, Method, dimensions, tElapsed, CPUtimeAtStart, cputime); 
fprintf( fInfo, 'End ground state sims: %s\n', datestr(endTime, 'dd mmm yy @ HH:MM:SS')); % end time

% save workspace to workspace folder
fnameWspace = strcat('xWorkspace-data/', fnameBase); % = 'GPEpaper#_Example#_#D'
save(strcat(fnameWspace, '-workspace-groundstate') )
save(strcat(fnameWspace, '-outputs-groundstate') , 'Phi_1', 'Outputs')

% (print) starting time for dynamical simulation
tStartDyn = tic; % starting time for dynamical sims timer
CPUtimeAtStartDyn = cputime; % CPU start time
fprintf( fInfo, 'Start dynamics: %s\n', datestr(now, 'dd mmm yy @ HH:MM:SS'));

%%
%{
DYNAMICAL SIMULATION
%}

%% Simulation methods

Computation = 'Dynamic';
Ncomponents = 1;
Type = 'Splitting';
Deltat = 1e-3;
Stop_time = 1;
Stop_crit = [];

Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);

%% We keep Geometry2D as-is

%% Physics2D: we add a random potential

X_0 = 0;
Y_0 = 0;
d = 4;
V_0 = 2;
Brownian = Brownian_Process2d(Method);
Physics2D = StochasticPotential_Var2d(Method, Physics2D, ...
    @(W,X,Y) V_0 * exp( -( (X-X_0).^2 + (Y-Y_0).^2 )/2*d^2 ).*W, [], ...
    @(t,X,Y) Brownian(t));

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

%% end
