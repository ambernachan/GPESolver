%{
Low-interactions system without rotations in 3D
%}

function [] = trial_Gaussian3D_weakinteractions_Delta(info)

    close all;

    %% Setting variables
    % Determine interaction strength compared to kinetic energy
    S = info.params.S;
    % Setting simulation space
    xlim = info.params.boxlimits(1); ylim = info.params.boxlimits(2); zlim = info.params.boxlimits(3);
    %Nx = Ngridpts(1); Ny = Ngridpts(2); Nz = Ngridpts(3);
    Nx = info.params.Ngridpts; Ny = info.params.Ngridpts; Nz = info.params.Ngridpts;
    % Setting physical parameters
    Delta = info.params.delta;
    gx = info.params.gammas(1); gy = info.params.gammas(2); gz = info.params.gammas(3);
    
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
    Stop_crit = {'MaxNorm', gx*1e-10};
    Max_iter = 1e3;

    Method = Method_Var3d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit, Max_iter);

    %% Geometry3D
    
    xmin = -xlim;
    xmax = xlim;
    ymin = -ylim;
    ymax = ylim;
    zmin = -zlim;
    zmax = zlim;
    
    Geometry3D = Geometry3D_Var3d(xmin, xmax, ymin, ymax, zmin, zmax, Nx, Ny, Nz);

    %% Take into account possible higher # of chis/S
    
    %info = Info(name_from_filename(mfilename), dimensions);
    
    % Save ground state method
    Method_ground = Method;
    save(info.get_workspace_path('fittingdata'),'Method_ground')
    clear Method_ground;
    
    

    %% Physics3D

    Beta = 4*pi*S;
    Omega = 0;
    Physics3D = Physics3D_Var3d(Method, Delta, Beta, Omega);
    Physics3D = Dispersion_Var3d(Method, Physics3D); % !!!
    %Physics3D = Potential_Var3d(Method, Physics3D); % std quadratic potential
    %Physics3D = Potential_Var3d(Method, Physics3D, @(X,Y,Z) quadratic_potential3d(1, 5, 1, X, Y, Z));

    Physics3D = Potential_Var3d(Method, Physics3D, @(X,Y,Z) quadratic_potential3d(gx, gy, gz, X, Y, Z));
    Physics3D = Nonlinearity_Var3d(Method, Physics3D); % std cubic nonlinearity
    %Physics3D = Gradientx_Var3d(Method, Physics3D, @(x,y) -1i*Omega*y);
    %Physics3D = Gradienty_Var3d(Method, Physics3D, @(x,y) 1i*Omega*x);

    %% Defining a starting function Phi_0

    InitialData_choice = 1; % Gaussian initial data

    % find w3d
    exp = ExpPhis(S, [gx,gy,gz], [0,0,0], Geometry3D);
    W = exp.getW();
    w = W.d3;
    clear exp W

    sigma_w_delta = struct();
    sigma_w_delta.x = w * (Delta/2)^(0.25) / sqrt(gx); sigma_w_delta.x = sigma_w_delta.x(sigma_w_delta.x>0);
    sigma_w_delta.y = w * (Delta/2)^(0.25) / sqrt(gy); sigma_w_delta.y = sigma_w_delta.y(sigma_w_delta.y>0);
    sigma_w_delta.z = w * (Delta/2)^(0.25) / sqrt(gz); sigma_w_delta.z = sigma_w_delta.z(sigma_w_delta.z>0);

    X0 = 0;
    Y0 = 0;
    Z0 = 0;
    gamma_x = 1;
    gamma_y = 1;
    gamma_z = 1;

    Phi_0 = InitialData_Var3d(Method, Geometry3D, Physics3D, InitialData_choice, X0, Y0, Z0, gamma_x, gamma_y, gamma_z);

    %% version mgmt
    curdir = strsplit(pwd, '/');
    curdir = strsplit(curdir{end}, '\');
    curdir = curdir{end};

    %% Printing interaction strength

    info.add_info_separator();
    info.add_custom_info('Folder: \t %s \n', curdir); % print current folder
    info.add_info_separator();
    info.add_custom_info('S \t=\t %f \n', S); % print interaction strength
    info.add_custom_info('Beta \t=\t %f \n', Beta); % print interaction parameter Beta
    info.add_custom_info('w \t=\t %f \n', w); % print Gaussian parameter w
    info.add_custom_info('sigmas \t=\t [%.4f,%.4f,%.4f] \n', ...
        sigma_w_delta.x,sigma_w_delta.y,sigma_w_delta.z); % print Gaussian parameter sigma_xyz with Delta
    info.add_custom_info('Delta \t=\t %f \n', Delta); % print kinetic energy parameter Delta
    info.add_custom_info('gammas \t=\t [%.4f,%.4f,%.4f] \n', gx,gy,gz); % print gammas
    info.add_info_separator();

    %% Determining outputs

    Outputs = OutputsINI_Var3d(Method);

    %% Printing preliminary outputs

    Printing = 1;
    Evo = 15;
    Draw = 1;
    Print = Print_Var3d(Printing, Evo, Draw);

    %% RUN THE SIMULATION to find the ground state

    info.add_simulation_info(Geometry3D);
    [Phi_1, Outputs] = GPELab3d(Phi_0, Method, Geometry3D, Physics3D, Outputs, [], Print);

    %% Save the workspace & simulation info

    % save information about final iteration in info file
    info.add_result_info(Method, Outputs);
    % save workspace to workspace folder
    save(info.get_workspace_path('groundstate'));

    %% Draw & save solution

    close all;
    pause(2) % pauses the program for 2 seconds

    Draw_solution3d(Phi_0, Method, Geometry3D, Figure_Var3d());

    info.save_figure(1, 'initialdata', 'psi_sq',[],[]);
    info.save_figure(2, 'initialdata', 'angle',[],[]);

    Draw_solution3d(Phi_1, Method, Geometry3D, Figure_Var3d());

    info.save_figure(1, 'groundstate', 'psi_sq',[],[]);
    info.save_figure(2, 'groundstate', 'angle',[],[]);

    simulation_finished = 'Ground state simulation finished';

    save(info.get_workspace_path('fittingdata'),'simulation_finished', '-append')
    
    %%
    if info.params.dyn_simu
        %{
        DYNAMICAL SIMULATION
        %}

        %% Simulation methods

        Computation = 'Dynamic';
        Ncomponents = 1;
        Type = 'Splitting';
        Deltat = 1e-4;
        %Stop_time = 1;
        Max_iter = 200;
        %Stop_crit = [];
        Stop_crit = {'MaxNorm', 1e-12};

        %Method = Method_Var3d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit, Max_iter);
        Method = Method_Var3d(Computation, Ncomponents, Type, Deltat, [], Stop_crit, Max_iter);

        % Saving workspace with relevant data for fitting
        Method_dynamical = Method;
        save(info.get_workspace_path('fittingdata'), 'Method_dynamical', '-append');
        clear Method_dynamical;

        %% We keep Geometry3D as-is

        %% We keep Physics3D as-is

        %% Determining outputs

        Save_solution = 1;
        Outputs = OutputsINI_Var3d(Method, Save_solution);

        %% Printing preliminary outputs
        Printing = 1;
        Evo = 10;
        Draw = 0;
        Print = Print_Var3d(Printing, Evo, Draw);

        %% RUN THE SIMULATION to find the ground state

        [Phi, Outputs] = GPELab3d(Phi_1, Method, Geometry3D, Physics3D, Outputs, [], Print);

        %% Save the workspace & simulation info

        % save information about final simulation iteration in info file
        info.add_result_info(Method, Outputs);
        info.finish_info();
        save(info.get_workspace_path('dynamics'))

        %% Draw & save solution

        close all;
        pause(2) % pauses the program for 2 seconds

        Draw_solution3d(Phi, Method, Geometry3D, Figure_Var3d());

        info.save_figure(1, 'dynamics', 'psi_sq',[],[]);
        info.save_figure(2, 'dynamics', 'angle',[],[]);
        
        % Save dynamic sim results to workspace
        phi_dyn = PhiData(Phi, Geometry3D);
        save(info.get_workspace_path('fittingdata'), 'phi_dyn', '-append');
    end

    %% Save PhiData structures for fitting w/o the whole workspace

    phi_input = PhiData(Phi_0, Geometry3D);
    phi_ground = PhiData(Phi_1, Geometry3D);
   
    % Saving workspace with relevant data for fitting
    save(info.get_workspace_path('fittingdata'), 'gx', 'gy', 'gz', ...
        'phi_ground', 'phi_input', 'info', 'Method', ... % necessary data
        'S', 'w', 'Beta', 'Delta', 'sigma_w_delta', ... % additional data
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
    load(wspath)

    close all;
    %%%%%%%%%%%%%%%%%%%%% % set x boundaries
    if S >= 0.01
        limits = 3;
    elseif S < 0.01
        limits = [3, 1, 0.2, 0.1];
    end
    
    saveBunchOfPlots(wspath, limits)
    
    simulation_finished = 'yes!';
    save(wspath,'simulation_finished', '-append')

    %% end
end
