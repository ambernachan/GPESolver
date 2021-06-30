 %{
23Na in the F=1 manifold
Spinor BEC in 3D
No rotations, no magnetic field, quadratic optical trap
%}

function [] = spinor_GPE3D_dynamics(info, Phi_in)
    
    close all;
   
    %% Setting variables
    % Determine interaction strength compared to kinetic energy
%     atom = 'Na';
    if isfield(info.params, 'atom')
        atom = info.params.atom;
    else
        atom = 'Rb';
%         atom = 'Na';
        info.params.atom = atom;
        sprintf('Warning: atom set to default (87Rb) as atom type was not specified')
%         sprintf('Warning: atom set to default (23Na) as atom type was not specified')
    end
    if isfield(info.params, 'a0')
        a0 = info.params.a0;
    else
        a0 = getsimconst(['a0_' atom]);
        info.params.a0 = a0;
    end
    if isfield(info.params, 'a2')
        a2 = info.params.a2;
    else
        a2 = getsimconst(['a2_' atom]);
        info.params.a2 = a2;
    end
    N = getsimconst('N'); % number of particles
    hbar = getphysconst('hbar'); % in kg m^2 / s
    trapfreq = getsimconst('trap_freq'); % Trap strength in Hz (symmetric for now)
    atom_mass = getsimconst(['mass_' atom]); % Atom mass in kg
    spin_pair = getsimconst('spin_pair'); % hyperfine spin manifold (=1)
    info.params.spin_pair = spin_pair;
    
    aho = sqrt(hbar / (atom_mass * trapfreq));
    an = (2*a2+a0)/3;
    as = (a2-a0)/3;
    chin = N*an/aho;
    chis = N*as/aho;
    
    % Setting simulation space
    xlim = info.params.boxlimits(1); ylim = info.params.boxlimits(2); zlim = info.params.boxlimits(3);
    Nx = info.params.Ngridpts; Ny = info.params.Ngridpts; Nz = info.params.Ngridpts;
    
    % Setting physical parameters
    if isfield(info.params, 'delta')
        Delta = info.params.delta;
    else
        Delta = 0.5;
        info.params.delta = Delta;
    end
    if isfield(info.params, 'gammas')
        gx = info.params.gammas(1); gy = info.params.gammas(2); gz = info.params.gammas(3);
    else
        gx = 1; gy = 1; gz = 1;
        info.params.gammas(1) = gx;
        info.params.gammas(1) = gy;
        info.params.gammas(1) = gz;
    end
    if isfield(info.params, 'dimensions')
        dimensions = info.params.dimensions;
    else
        dimensions = 3;
        info.params.dimensions = dimensions;
    end
    
    % xi is the healing length; we derive a xi-n and a xi-s for self- and
    % spin-mixing interactions
    allthechis = [chin chis];
    alltheAs = [an as];
    XI = findhealinglengths(allthechis, alltheAs, atom);
    
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
%     Deltat = 0.1*dx^2;
    Deltat = 0.1*dx^2;
    Stop_time = 100;
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
    Beta = 1; % multiplication factor for Beta_n and Beta_s
    Betan = 4*pi*chin;
    Betas = 4*pi*chis;
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
    info.add_custom_info('Beta_n \t=\t %f \n', Betan); % print interaction parameter Beta_n
    info.add_custom_info('Beta_s \t=\t %f \n', Betas); % print interaction parameter Beta_s
    info.add_custom_info('Delta \t=\t %f \n', Delta); % print kinetic energy parameter Delta
    info.add_custom_info('gammas \t=\t [%.4f,%.4f,%.4f] \n', gx,gy,gz); % print gammas
    info.add_info_separator();
    
    %% Determining outputs

    Evo_outputs = 10; % Must be equal to or smaller than Evo from Print
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
        Directional_Magnetization(Method, Geometry3D, Phi, 'x', X, Y, Z, FFTX, FFTY, FFTZ);
    globaluserdef_outputs{6} = @(Phi,X,Y,Z,FFTX,FFTY,FFTZ) ...
        Directional_Magnetization(Method, Geometry3D, Phi, 'y', X, Y, Z, FFTX, FFTY, FFTZ);
    globaluserdef_outputs{7} = @(Phi,X,Y,Z,FFTX,FFTY,FFTZ) ...
        Directional_Magnetization(Method, Geometry3D, Phi, 'z', X, Y, Z, FFTX, FFTY, FFTZ);
    globaluserdef_names{1} = 'Longitudinal magnetization';
    globaluserdef_names{2} = 'Population psi+';
    globaluserdef_names{3} = 'Population psi0';
    globaluserdef_names{4} = 'Population psi-';
    globaluserdef_names{5} = 'Magnetization Mx';
    globaluserdef_names{6} = 'Magnetization My';
    globaluserdef_names{7} = 'Magnetization Mz';
    
    Outputs = OutputsINI_Var3d(Method, Evo_outputs, Save_solution, [], [], ...
        globaluserdef_outputs,globaluserdef_names);
    
    %% Printing preliminary outputs
    
    Printing = 1;
    Evo = 25;
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
    
    save(info.get_workspace_path('phi_ini'), 'Phi_in', 'Phi', 'Outputs', '-v7.3')

    %% Save the workspace & simulation info

    % save information about final simulation iteration in info file
    info.add_result_info(Method, Outputs);

    % saving dynamics workspace in v7.3 MAT file
    save(info.get_workspace_path('dynamics_v7.3'), '-v7.3');
    
    % saving dynamics workspace in 'normal' file
    save(info.get_workspace_path('dynamics'))
        
    %% Draw user-defined functions populations & magnetization
    
    close all;
    pause(2) % pauses the program for 2 seconds
    
    its = Outputs.Iterations;
    
    % Plot longitudinal magnetization
    plot_longmagnetization(its, Outputs.User_defined_global{1}, info, Outputs.Evo_outputs)
    % Plot population fractions
    plot_populationfractions(its, Outputs.User_defined_global(2:4), info, Outputs.Evo_outputs)
    % Plot population distribution on x-axis
    plot_populationdistribution(Geometry3D, Phi, info)
    % Plot magnetization distribution on x-axis
    plot_magnetizationdistribution(Geometry3D, Phi, info)
    % Plot transverse & longitudinal magnetization
%     plot_magnetizations(its, Outputs.User_defined_global{5}, info, Outputs.Evo_outputs)
    
    %% Time plots
    % magnetization distribution on x-axis
    timeslider_magnetizationdistribution(Geometry3D, Outputs.Solution, info)
    % population distribution on x-axis
    timeslider_populationdistribution(Geometry3D, Outputs.Solution, info)
    
    
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
    timeunits = 1/trapfreq;
   
    % Saving workspace with relevant data for fitting
    save(info.get_workspace_path('fittingdata'), 'gx', 'gy', 'gz', ...
        'phi_ground', 'phi_input', 'info', 'Method', ... % necessary data
        'Beta', 'Betan', 'Betas', 'Delta', 'Ncomponents', ... % additional data
        'chin', 'chis', 'trapfreq', 'timeunits', ... % additional data
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
    
    limits = [5,3,1];
    
    saveBunchOfPlots(wspath, limits)
    
    simulation_finished = 'yes!';
    save(wspath,'simulation_finished', '-append')

    info.finish_info();
    
    %% end
end
