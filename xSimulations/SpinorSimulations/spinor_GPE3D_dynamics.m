 %{
23Na in the F=1 manifold
Spinor BEC in 3D
No rotations, no magnetic field, quadratic optical trap
%}

function [] = spinor_GPE3D_dynamics(info, params)
    
    close all;
    Phi_in = params.Phi_input;
   
    %% Setting variables
    
    % Setting simulation space
    xlim = info.params.boxlimits(1); ylim = info.params.boxlimits(2); zlim = info.params.boxlimits(3);
    Nx = info.params.Ngridpts; Ny = info.params.Ngridpts; Nz = info.params.Ngridpts;
    
    % Setting physical parameters
    if isprop(info.params, 'delta')
        Delta = info.params.delta;
    else
        Delta = 0.5;
        info.params.delta = Delta;
    end
    if isprop(info.params, 'gammas')
        gx = info.params.gammas(1); gy = info.params.gammas(2); gz = info.params.gammas(3);
    else
        gx = 1; gy = 1; gz = 1;
        info.params.gammas(1) = gx;
        info.params.gammas(1) = gy;
        info.params.gammas(1) = gz;
    end
    if isprop(info.params, 'dimensions')
        dimensions = info.params.dimensions;
    else
        dimensions = 3;
        info.params.dimensions = dimensions;
    end
    
    % xi is the healing length; we derive a xi-n and a xi-s for self- and
    % spin-mixing interactions
    allthechis = [info.params.chin info.params.chis];
    alltheAs = [info.params.an info.params.as];
    XI = findhealinglengths(allthechis, alltheAs, info.params.atom);
    
    %% Run a ground state simulation if Phi_in is not yet given.
    
    if ~exist('Phi_in', 'var')
        groundst_simulation = 'spinor_GPE3D_ground';
        ground_info = multi_runner(groundst_simulation, info.params.boxlimits, info.params.Ngridpts);
        ws_ground = ground_info.get_workspace_path('groundstate'); % ground state workspace path/name
        ws_ground = load(ws_ground); % ground state workspace struct of variables
        names = fieldnames(ws_ground); % variable names cell in ground state ws
        for i = 1:length(names)
            if strcmp(names{i}, 'Phi_1')
                Phi_in = ws_ground.(names{21});
            end
        end
        if ~exist('Phi_in', 'var')
            error('Something went wrong with the ground state simulation, please check.')
        end
    end
        
    %{
    DYNAMICAL SIMULATION
    %}

    %% Simulation methods

    Computation = 'Dynamic';
    Ncomponents = 3;
    Type = 'Splitting'; % defaults to Strang splitting
%     Type = 'Relaxation';
    dx = (2*xlim / (Nx-1));
    Deltat = info.params.dt;
%     Deltat = 0.1*(dx)^2;
    LimitingIter = 120000; % A limitation to the iterations bc time
    Stop_time = floor(min(100, (LimitingIter*Deltat)));
    Stop_crit = [];
    Max_iter = 10000;
%         Stop_crit = {'MaxNorm', 1e-12};
%         Precond_type = 'FLaplace'; % defaults to 'FLaplace'
%         Precond_type = []; % defaults to 'FLaplace'
%     Precond_type = 'FThomasFermi'; % defaults to 'FLaplace'
%     Output = 1;

    Method = Method_Var3d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit, Max_iter);
%     Method = Method_Var3d(Computation, Ncomponents, Type, Deltat, Stop_time, [], [], Precond_type, Output, Splitting);

    % Saving workspace with relevant data for fitting
    Method_dynamical = Method;
    
    %% Save dynamical method to fittingdata workspace
    if isfile(info.get_workspace_path('fittingdata'))
        save(info.get_workspace_path('fittingdata'), 'Method_dynamical', '-append');
    else
        save(info.get_workspace_path('fittingdata'), 'Method_dynamical');
    end
    clear Method_dynamical;

    %% Geometry3D
    
    xmin = -xlim;   xmax = xlim;
    ymin = -ylim;   ymax = ylim;
    zmin = -zlim;   zmax = zlim;
    
    Geometry3D = Geometry3D_Var3d(xmin, xmax, ymin, ymax, zmin, zmax, Nx, Ny, Nz);%% Geometry3D
        
    %% Physics3D

    % Delta is already defined
    Beta = info.params.beta; % multiplication factor for Beta_n and Beta_s
    Betan = info.params.betan;
    Betas = info.params.betas;
    Omega = 0;
    Physics3D = Physics3D_Var3d(Method, Delta, Beta, Omega);
    Physics3D = Dispersion_Var3d(Method, Physics3D); % !!!
    Physics3D = Potential_Var3d(Method, Physics3D); % std quadratic potential
%     Physics3D = Potential_Var3d(Method, Physics3D, @(X,Y,Z) quadratic_potential3d(1, 5, 1, X, Y, Z));
%     Physics3D = Potential_Var3d(Method, Physics3D, @(X,Y,Z) quadratic_potential3d(gx, gy, gz, X, Y, Z));
    
%     Physics3D = Nonlinearity_Var3d(Method, Physics3D, Coupled_Cubic3d_spin1(Betan,Betas)); % cubic nonlinearity with off-diagonal coupling
    Physics3D = Nonlinearity_Var3d(Method, Physics3D, Coupled_Cubic3d_spin1(Betan,Betas), [], ...
        Coupled_CubicEnergy3d_spin1(Betan,Betas)); % cubic nonlinearity with off-diagonal coupling
    
    %Physics3D = Gradientx_Var3d(Method, Physics3D, @(x,y) -1i*Omega*y);
    %Physics3D = Gradienty_Var3d(Method, Physics3D, @(x,y) 1i*Omega*x);
    
    %% version mgmt
    curdir = strsplit(pwd, '/');
    curdir = strsplit(curdir{end}, '\');
    curdir = curdir{end};

    %% Printing interaction strength

    info.add_info_separator();
    info.add_custom_info('Folder: \t %s \n', curdir); % print current folder
    info.add_info_separator();
    info.add_custom_info('atom \t=\t %s \n', info.params.atom);
    info.add_custom_info('Beta_n \t=\t %f \n', info.params.betan); % print interaction parameter Beta_n
    info.add_custom_info('Beta_s \t=\t %f \n', info.params.betas); % print interaction parameter Beta_s
%     info.add_custom_info('Delta \t=\t %f \n', Delta); % print kinetic energy parameter Delta
%     info.add_custom_info('gammas \t=\t [%.4f,%.4f,%.4f] \n', gx,gy,gz); % print gammas
    info.add_custom_info('dt \t=\t %f \n', info.params.dt);
    info.add_custom_info('dx \t=\t %f \n', dx);
    info.add_info_separator();
    
    %% Determining outputs

    % Must be equal to or smaller than Evo from Print
    Evolim = round((3*(Nx)^3*Stop_time / (info.params.dt*7e7))/5)*5;
    Evo_outputs = max(10, Evolim);
    Save_solution = 1;
    
    globaluserdef_outputs{1} = @(Phi,X,Y,Z,FFTX,FFTY,FFTZ) ...
        Magnetization(Method, Geometry3D, Phi, X, Y, Z, FFTX, FFTY, FFTZ);
    globaluserdef_outputs{2} = @(Phi,X,Y,Z,FFTX,FFTY,FFTZ) ...
        Population(Method, Geometry3D, Phi, X, Y, Z, FFTX, FFTY, FFTZ, 1);
    globaluserdef_outputs{3} = @(Phi,X,Y,Z,FFTX,FFTY,FFTZ) ...
        Population(Method, Geometry3D, Phi, X, Y, Z, FFTX, FFTY, FFTZ, 0);
    globaluserdef_outputs{4} = @(Phi,X,Y,Z,FFTX,FFTY,FFTZ) ...
        Population(Method, Geometry3D, Phi, X, Y, Z, FFTX, FFTY, FFTZ, -1);
    globaluserdef_outputs{5} = @(Phi,X,Y,Z,FFTX,FFTY,FFTZ) ...
        Directional_Magnetization(Method, Geometry3D, Phi, 'Mx', X, Y, Z, FFTX, FFTY, FFTZ);
    globaluserdef_outputs{6} = @(Phi,X,Y,Z,FFTX,FFTY,FFTZ) ...
        Directional_Magnetization(Method, Geometry3D, Phi, 'My', X, Y, Z, FFTX, FFTY, FFTZ);
    globaluserdef_outputs{7} = @(Phi,X,Y,Z,FFTX,FFTY,FFTZ) ...
        Directional_Magnetization(Method, Geometry3D, Phi, 'Mz', X, Y, Z, FFTX, FFTY, FFTZ);
    globaluserdef_outputs{8} = @(Phi,X,Y,Z,FFTX,FFTY,FFTZ) ...
        Directional_Magnetization(Method, Geometry3D, Phi, 'x', X, Y, Z, FFTX, FFTY, FFTZ);
    globaluserdef_outputs{9} = @(Phi,X,Y,Z,FFTX,FFTY,FFTZ) ...
        Directional_Magnetization(Method, Geometry3D, Phi, 'y', X, Y, Z, FFTX, FFTY, FFTZ);
    globaluserdef_outputs{10} = @(Phi,X,Y,Z,FFTX,FFTY,FFTZ) ...
        Directional_Magnetization(Method, Geometry3D, Phi, 'z', X, Y, Z, FFTX, FFTY, FFTZ);
    globaluserdef_names{1} = 'Longitudinal magnetization';
    globaluserdef_names{2} = 'Population psi+';
    globaluserdef_names{3} = 'Population psi0';
    globaluserdef_names{4} = 'Population psi-';
    globaluserdef_names{5} = 'Magnetization Mx';
    globaluserdef_names{6} = 'Magnetization My';
    globaluserdef_names{7} = 'Magnetization Mz';
    globaluserdef_names{8} = 'Magnetization Fx';
    globaluserdef_names{9} = 'Magnetization Fy';
    globaluserdef_names{10} = 'Magnetization Fz';
    
    Outputs = OutputsINI_Var3d(Method, Evo_outputs, Save_solution, [], [], ...
        globaluserdef_outputs,globaluserdef_names);
    
    %% Printing preliminary outputs
    
    Printing = 1;
    % Must be equal to or bigger than Evo_outputs
    Evo = max(25, Evo_outputs);
    Draw = 0;
    Print = Print_Var3d(Printing, Evo, Draw);

    %% RUN THE SIMULATION to find the ground state

%         A = [0.1 0.2 0.3];
%         d = [0.6 0.5 0.4];
%         for i = 1:3
%             Random_Phase{i} = Stationary_Gaussian_Field3d(Geometry3D,...
%             @(X,Y,Z) A(i)*exp(-(X.^2 + Y.^2 + Z.^2)/(2*d(i)^2)));
%             testerPhi{i} = exp(-2i*pi*Random_Phase{i});
%         end
%         [Phi, Outputs] = GPELab3d(testerPhi, Method, Geometry3D, Physics3D, Outputs, [], Print);

%         Phi_i{3} = Phi_1{3} ./ sqrt(2);
%         Phi_i{2} = Phi_1{3} ./ sqrt(2);
%         Phi_i{1} = Phi_1{1};

    info.add_simulation_info(Geometry3D);
%         [Phi, Outputs] = GPELab3d(Phi_i, Method, Geometry3D, Physics3D, Outputs, [], Print);
    [Phi, Outputs] = GPELab3d(Phi_in, Method, Geometry3D, Physics3D, Outputs, [], Print);
    
    save(info.get_workspace_path('phi_ini'), 'Phi_in', 'Phi', 'Outputs', 'info', '-v7.3')
    
    %% Save the workspace & simulation info

    % save information about final simulation iteration in info file
    info.add_result_info(Method, Outputs);

    % saving dynamics workspace in v7.3 MAT file
    save(info.get_workspace_path('dynamics_v7.3'), '-v7.3');
    
    % saving dynamics workspace in 'normal' file
    save(info.get_workspace_path('dynamics'))
    
    %% create file that shows type of atom and simulation stats
    fname = createTextFileName(info, Geometry3D, Method, Outputs.Iterations*Outputs.Evo_outputs);
    input_str = 'XXX';
    if exist('inputFunctionDate', 'var')
        input_str = inputFunctionDate;
    end
    fname = [fname '_fromstateXXX'  '_usinginputf' input_str];
    fname = [info.fulldir '/' fname '.txt'];
    fileID = fopen(fname,'w');
    fclose(fileID);
    
    %% Draw user-defined functions populations & magnetization
    
    close all;
    pause(2) % pauses the program for 2 seconds
    
    its = Outputs.Iterations;
    
    M = Outputs.User_defined_global(5:7);
    F = Outputs.User_defined_global(8:10);
    
    % Plot longitudinal magnetization
    plot_longmagnetization(its, Outputs.User_defined_global{1}, info, Outputs.Evo_outputs)
    % Plot population fractions
    plot_populationfractions(its, Outputs.User_defined_global(2:4), info, Outputs.Evo_outputs, Method)
    % Plot population distribution on x-axis
    plot_populationdistribution(Geometry3D, Phi, info)
    % Plot magnetization distribution on x-axis
    plot_magnetizationdistribution(Geometry3D, Phi, info)
    % Plot transverse & longitudinal magnetization
    plot_magnetizations(its, F, info, Outputs.Evo_outputs, Method)
    
    %% Time plots
    % magnetization distribution on x-axis
    timeslider_magnetizationdistribution(Geometry3D, Outputs.Solution, info)
    % population distribution on x-axis
    timeslider_populationdistribution(Geometry3D, Outputs.Solution, info)
    % phase distribution as a sliced 3d function
    timeslider_phase(Geometry3D, Outputs.Solution, info)
    
    % adding extra parameters to workspace
    save(info.get_workspace_path('dynamics_v7.3'), 'M', 'F', 'its', '-append');
    
    %% Draw & save solution

    close all;
    pause(2) % pauses the program for 2 seconds

    % Set figure names
    for i = 1 : Ncomponents
        phistr = 'phi_sq';
        anglestr = 'angle';
        if Ncomponents > 1
            component_str = componentstr(Ncomponents); % ['+0-']
            phisqname{i} = [phistr '_' component_str(i)];
            anglename{i} = [anglestr '_' component_str(i)];
            if i == Ncomponents
                phisqname{i+1} = [phistr '_tot'];
            end
        else
            phisqname = {phistr};
            anglename = {anglestr};
        end
    end
    
    Draw_solution3d(Phi, Method, Geometry3D, Figure_Var3d());
    for i = 1 : Ncomponents
        info.save_figure(2*i-1, 'dynamics', phisqname{i},[],[]);
        info.save_figure(2*i, 'dynamics', anglename{i},[],[]);

        if ((Ncomponents > 1) && (i == Ncomponents))
            info.save_figure(2*i+1, 'dynamics', phisqname{i+1},[],[]);
        end
    end

    % Save dynamic sim results to workspace
    phi_dyn = PhiData(Phi, Geometry3D);
    save(info.get_workspace_path('fittingdata'), 'phi_dyn', '-append');

    %% Save PhiData structures for fitting w/o the whole workspace

%     phi_initial = PhiData(Phi_0, Geometry3D);
    phi_input = PhiData(Phi_in, Geometry3D);
    phi_ground = PhiData(Phi_in, Geometry3D);
    timeunits = 1/info.params.trapfreq;
   
    % Saving workspace with relevant data for fitting
    save(info.get_workspace_path('fittingdata'), 'gx', 'gy', 'gz', ...
        'phi_ground', 'phi_input', 'info', 'Method', 'timeunits',... % necessary data
        '-append'); % to not overwrite Method_ground

    simulation_finished = 'yes, but no plots yet!';
    save(info.get_workspace_path('fittingdata'),'simulation_finished', '-append')
    
    %% Some new output

    fulldir = info.fulldir;
    wspath = info.get_workspace_path('fittingdata');
    datestring = info.creationTimeString;
    save(info.get_workspace_path('fittingdata'),'fulldir', 'wspath', 'datestring', 'info', '-append')
    
    clearvars -except wspath
    analysisPhi3d(wspath);

    close all;
    
    limits = [8];
    saveBunchOfPlots(wspath, limits)
    
    simulation_finished = 'yes!';
    save(wspath,'simulation_finished', '-append')
    
    load(wspath, 'info')
    info.finish_info();
    
    %% end
end
