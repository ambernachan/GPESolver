%{
Collision of two bright solitons for a system of GPEs in 1d
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

Computation = 'Dynamic';
Ncomponents = 2;
Type = 'Relaxation';
Deltat = 1e-2;
Stop_time = 25;

Method = Method_Var1d(Computation, Ncomponents, Type, Deltat, Stop_time);

%% Geometry1D

xmin = -40;
xmax = 40;
Nx = 2^11 + 1;
Geometry1D = Geometry1D_Var1d(xmin, xmax, Nx);

%% Physics1D with coupled nonlinearity and default Laplacian operator

Delta = 1;
Beta = 1;
alpha_1 = 0.25;
alpha_2 = -0.1965;
Physics1D = Physics1D_Var1d(Method, Delta, Beta);
% using the default Laplacian dispersion operator:
Physics1D = Dispersion_Var1d(Method, Physics1D);

Coupled_NL{1,1} = @(Phi, X) alpha_1*abs(Phi{1}).^2 + (alpha_1+2*alpha_2)*abs(Phi{2}).^2;
Coupled_NL{2,1} = @(Phi, X) 0;
Coupled_NL{1,2} = @(Phi, X) 0;
Coupled_NL{2,2} = @(Phi, X) alpha_1*abs(Phi{2}).^2 + (alpha_1+2*alpha_2)*abs(Phi{1}).^2;

Physics1D = Nonlinearity_Var1d(Method, Physics1D, Coupled_NL);

%% Defining a starting function Phi_0

c_1 = 1.2;
n_1 = 0.03;
b_1 = sqrt(n_1+(c_1^2)/4);
X_1 = -15;
X = Geometry1D.X;

Phi_0{1} = sqrt(2/abs(alpha_1)) * b_1 * sech(b_1*(X-X_1)).* exp(1i*c_1*X/2+n_1);
c_r = -0.5;
n_r = 0.1;
b_r = sqrt(n_r + (c_r^2)/4);
X_r = 0;
X = Geometry1D.X;
Phi_0{2} = sqrt(2/abs(alpha_1)) * b_r * sech(b_r*(X-X_r)).* exp(1i*c_r*X/2+n_r);

%% Save the position of the soliton in a user-defined output function
Solution_save = 1;
Outputs_iterations = 10;
Output_function{1} = @(Phi, X, FFTX) Geometry1D.dx * sum(X.*abs(Phi).^2);
Output_name{1} = 'Position of the soliton';
Outputs = OutputsINI_Var1d(Method, Outputs_iterations, Solution_save, Output_function, Output_name);

%% Printing preliminary outputs
Printing = 1;
Evo = 15;
Draw = 1;
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

savefig(1, strcat(fnameFigspace, '-psi_sq_sq-1'))
savefig(2, strcat(fnameFigspace, '-psi_sq_sq-2'))
savefig(3, strcat(fnameFigspace, '-angle-1'))
savefig(4, strcat(fnameFigspace, '-angle-2'))

Draw_Timesolution1d(Outputs, Method, Geometry1D, Figure_Var1d);

savefig(1, strcat(fnameFigspace, '-psi_sq-1'))
savefig(2, strcat(fnameFigspace, '-psi_sq-2'))

%% Draw user-defined output function : position of the soliton

Time = [0 : 0.1 : 25];

figure(3)
plot(Time, Outputs.User_defined_local{1,1});
xlabel('Time')
ylabel('Position of $\Psi_1$')
savefig(strcat(fnameFigspace, '-position-psi1'))

figure(4)
plot(Time, Outputs.User_defined_local{1,1});
xlabel('Time')
ylabel('Position of $\Psi_2$')
savefig(strcat(fnameFigspace, '-position-psi2'))

%% end
