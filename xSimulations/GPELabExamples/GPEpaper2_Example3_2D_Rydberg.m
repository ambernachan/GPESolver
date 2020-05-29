%{
Dynamics of a Rydberg-dressed Bose-Einstein condensate in 2d
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
Deltat = 1e-2;
Stop_time = [];
Stop_crit = {'Energy', 1e-6};

Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);

%% Geometry2D

xmin = -15;
xmax = 15;
ymin = -15;
ymax = 15;
Nx = 2^9 + 1;
Ny = 2^9 + 1;
Geometry2D = Geometry2D_Var2d(xmin, xmax, ymin, ymax, Nx, Ny);

%% Physics2D for BEC with cubic+ potential and Rydberg nonlinear interaction

Delta = 0.5;
Beta = 20000;
Physics2D = Physics2D_Var2d(Method, Delta, Beta);

omega = 5;
V_0 = 750;
d = 0.7;
Physics2D = Potential_Var2d(Method, Physics2D, @(x,y) (omega^2/2)*(x.^2+y.^2) + ...
    V_0*exp(-x.^2/(2*d^2)));

R0 = 0.8;
Physics2D = Nonlinearity_Var2d(Method, Physics2D, @(phi,fftx,ffty) RydBergInteraction2d(R0,phi,fftx,ffty));

%% Defining a starting function Phi_0

InitialData_choice = 1; % initialization by a Gaussian function
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
fnameWspace = strcat('xWorkspace-data/', fnameBase);
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
Type = 'Relaxation';
Deltat = 1e-3;
Stop_time = 5;
Stop_crit = [];

Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);

%% We keep Geometry2D as-is

%% Physics2D: we remove the double-well potential

Delta = 0.5;
Beta = 20000;
Physics2D = Physics2D_Var2d(Method, Delta, Beta);

omega = 5;
Physics2D = Potential_Var2d(Method, Physics2D, @(x,y) (omega^2/2)*(x.^2+y.^2));

R0 = 0.8;
Physics2D = Nonlinearity_Var2d(Method, Physics2D, @(phi,fftx,ffty) RydBergInteraction2d(R0,phi,fftx,ffty));%, ...
    %[], @(phi, x, y, fftx, ffty) RydBergInteraction_energy2d(R0,phi,fftx, ffty));

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
