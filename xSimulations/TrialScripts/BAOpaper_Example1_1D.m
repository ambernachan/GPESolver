%{
Collision of two bright solitons for a system of GPEs in 1d
%}

clear;

%% Saving info files

workingFile = mfilename;
fnameBase = workingFile(1:20); %'BAOpaper_Example#_#D'%

fnameInfo = strcat('xOutputs/', fnameBase, '-info.txt'); % infofile name
fInfo = fopen(fnameInfo, 'w'); % create new file with simulation info
tStart = tic; % starting time for timer
CPUtimeAtStart = cputime; % CPU start time

fprintf( fInfo, 'Start: %s\n', datestr(now, 'dd mmm yy @ HH:MM:SS')); % print start time

%% Simulation methods

Computation = 'Dynamic';
Ncomponents = 1;
Type = 'Splitting';
Deltat = 1e-5;
Stop_time = 40;
%Stop_crit = {'MaxNorm', 5e-7};

Method = Method_Var1d(Computation, Ncomponents, Type, Deltat, Stop_time);

%% Geometry1D

xmin = -16;
xmax = 16;
Nx = 2^8 + 1;
Geometry1D = Geometry1D_Var1d(xmin, xmax, Nx);

%% Physics1D with coupled nonlinearity and default Laplacian operator

epsilon = 0.1;
kappa = 1.2649;

Delta = (1/2)*epsilon;
Beta = kappa/epsilon;
gamma_x = sqrt(1/epsilon);

Physics1D = Physics1D_Var1d(Method, Delta, Beta);
Physics1D = Dispersion_Var1d(Method, Physics1D);
Physics1D = Potential_Var1d(Method, Physics1D, @(X) quadratic_potential1d(gamma_x,X));
Physics1D = Nonlinearity_Var1d(Method, Physics1D, Cubic_w_Beta1d(Method, Beta));

%% Defining a starting function Phi_0

X = Geometry1D.X;

Phi_0{1} = 1/(pi*epsilon)^(1/4) * exp(-X.^2/(2*epsilon));

%% Save the position of the soliton in a user-defined output function
Solution_save = 1;
Outputs_iterations = 1000;
Output_function{1} = @(Phi, X, FFTX) Geometry1D.dx * sum(X.*abs(Phi).^2);
Output_name{1} = 'Position of the soliton';
Outputs = OutputsINI_Var1d(Method, Outputs_iterations, Solution_save, Output_function, Output_name);

%% Printing preliminary outputs
Printing = 1;
Evo = 1000;
Draw = 0;
Print = Print_Var1d(Printing, Evo, Draw);

%% RUN THE SIMULATION

[Phi, Outputs] = GPELab1d(Phi_0, Method, Geometry1D, Physics1D, Outputs, [], Print);

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

%% Draw & save solution
fnameFigspace = strcat('xFigures/', fnameBase);

%savefig(1, strcat(fnameFigspace, '-psi_sq_sq-1'))
%savefig(2, strcat(fnameFigspace, '-psi_sq_sq-2'))
%savefig(3, strcat(fnameFigspace, '-angle-1'))
%savefig(4, strcat(fnameFigspace, '-angle-2'))

Draw_Timesolution1d(Outputs, Method, Geometry1D, Figure_Var1d);

savefig(1, strcat(fnameFigspace, '-psi_sq-1'))

%% Draw user-defined output function : position of the soliton

Time = [0 : 1 : (Outputs.Iterations-1)];
%Time = [0 : 1 : 2500];
%Time = [0 : 0.1 : 25];

figure(3)
plot(Time, Outputs.User_defined_local{1,1});
xlabel('Time')
ylabel('Position of \Psi_1')
savefig(strcat(fnameFigspace, '-position-psi1'))

figure(4)
plot(Time, Outputs.x_rms{1});
xlabel('Time')
ylabel('\sigma_x (x_{RMS})')
savefig(strcat(fnameFigspace, '-xrms'))

%% end
