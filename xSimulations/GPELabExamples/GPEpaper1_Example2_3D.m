%{
Ground state calculation of a 3d GPE with a quadratic potential, a cubic
nonlinearity and a rotational operator
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

%% Method

Computation = 'Ground';
Ncomponents = 1;
Type = 'BESP';
Deltat = 5e-1;
Stop_time = [];
Stop_crit = {'MaxNorm', 1e-6};

Method = Method_Var3d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);

%% Geometry3D

xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;
zmin = -15;
zmax = 15;

Nx = 2^7 + 1;
Ny = 2^7 + 1;
Nz = 2^7 + 1;

Geometry3D = Geometry3D_Var3d(xmin,xmax, ymin,ymax, zmin,zmax, Nx,Ny,Nz);

%% Physics3D

Delta = 0.5;
Beta = 500;
Omega = [0,0,0.7];
gammax = 1;
gammay = 1;
gammaz = 0.5;

Physics3D = Physics3D_Var3d(Method,Delta,Beta,Omega);
Physics3D = Dispersion_Var3d(Method, Physics3D); % default Laplacian operator

Physics3D = Potential_Var3d(Method, Physics3D, @(X,Y,Z) ...
    quadratic_potential3d(gammax,gammay,gammaz, X,Y,Z));

Physics3D = Gradientx_Var3d(Method,Physics3D);
Physics3D = Gradienty_Var3d(Method,Physics3D);
Physics3D = Gradientz_Var3d(Method,Physics3D);

Physics3D = Nonlinearity_Var3d(Method,Physics3D);

%% Defining a starting function Phi_0

InitialData_choice = 2; % Thomas-Fermi approximation
Phi_0 = InitialData_Var3d(Method,Geometry3D,Physics3D,InitialData_choice);

%% Determining outputs

Outputs = OutputsINI_Var3d(Method);

%% Printing preliminary outputs

Printing = 1;
Evo = 15;
Draw = 1;
Print = Print_Var3d(Printing,Evo,Draw);

%% RUN THE SIMULATION

[Phi,Outputs] = GPELab3d(Phi_0,Method,Geometry3D,Physics3D,Outputs,[],Print);

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
