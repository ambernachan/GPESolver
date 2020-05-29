%{
Ground state calculation of a system of 2d GPEs modeling a spin-orbit-coupled
BEC under rotation
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

%% Simulation methods

Computation = 'Ground';
Ncomponents = 2;
Type = 'BESP';
Deltat = 5e-1;
Stop_time = [];
Stop_crit = {'MaxNorm', 1e-5};

Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);

%% Geometry2D

xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;
Nx = 2^9 + 1;
Ny = 2^9 + 1;

Geometry2D = Geometry2D_Var2d(xmin,xmax, ymin,ymax, Nx,Ny);

%% Physics2D

Delta = 0.5;
Beta = 1000;
Beta_coupled = [1,2; 2,1];
Omega = 0.1;
Kappa = 1.75;

Physics2D = Physics2D_Var2d(Method, Delta, Beta, Omega);
Physics2D = Dispersion_Var2d(Method, Physics2D, Dispersion_SpinOrbit2d(Kappa));
Physics2D = Potential_Var2d(Method, Physics2D);
Physics2D = Nonlinearity_Var2d(Method, Physics2D, ...
    Coupled_Cubic2d(Method, Beta_coupled), [], ...
    Coupled_Cubic_energy2d(Method, Beta_coupled));
Physics2D = Gradientx_Var2d(Method, Physics2D);
Physics2D = Gradienty_Var2d(Method, Physics2D);

%% Defining a starting function Phi_0

InitialData_choice = 2; % Thomas-Fermi approximation
Phi_0 = InitialData_Var2d(Method, Geometry2D, Physics2D, InitialData_choice);

%% Determining outputs

Outputs = OutputsINI_Var2d(Method);

%% Printing preliminary outputs

Printing = 1;
Evo = 15;
Draw = 1;
Print = Print_Var2d(Printing,Evo,Draw);

%% RUN THE SIMULATION

[Phi,Outputs] = GPELab2d(Phi_0,Method,Geometry2D,Physics2D,Outputs,[],Print);

%% Save the workspace & simulation info
endTime = now;
tElapsed = toc(tStart); % save time elapsed
dimensions = str2double(fnameBase(20)); % get dimensionality from file name

% save information about final simulation iteration in info file
add_info( fInfo, Outputs, Method, dimensions, tElapsed, CPUtimeAtStart, cputime); 
fprintf( fInfo, 'End: %s\n', datestr(endTime, 'dd mmm yy @ HH:MM:SS')); % end time
fclose( fInfo );

% save workspace to workspace folder
fnameWspace = strcat('xWorkspace-data/', fnameBase); % = 'GPEpaper#_Example#_#D'
save(strcat(fnameWspace, '-workspace') )
save(strcat(fnameWspace, '-outputs') , 'Phi', 'Outputs')

%% end
