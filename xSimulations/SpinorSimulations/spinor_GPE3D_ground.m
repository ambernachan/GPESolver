%{
23Na in the F=1 manifold
Spinor BEC in 3D
No rotations, no magnetic field, quadratic optical trap
One can change the type of atom by changing :
    - the atom mass
    - the a0 and a2 parameters
%}

function [] = spinor_GPE3D_ground(info)
    
    close all;
   
    %% Setting variables
    % Determine interaction strength compared to kinetic energy
%     atom = 'Na';
    if isfield(info.params, 'atom')
        atom = info.params.atom;
    else
        atom = 'Rb';
        info.params.atom = atom;
        sprintf('Warning: atom set to default (87Rb) as atom type was not specified')
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
    
    %%
    %{
    GROUND STATE SIMULATION TO GET STARTING FUNCTION PHI_0 FOR DYNAMIC SIMULATION
    %}

    %% Simulation methods

    Computation = 'Ground';
    Ncomponents = 3;
    Type = 'BESP';
    Deltat = 0.0015625;
    Stop_time = [];
    Stop_crit = {'MaxNorm', 1e-6};
    Max_iter = 2000;
%     Max_iter = 2500;

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
    
%     info = Info(name_from_filename(mfilename), dimensions);
    
%     Save ground state method
    Method_ground = Method;
    save(info.get_workspace_path('fittingdata'),'Method_ground')
    clear Method_ground;
    
    
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

    %% Defining a starting function Phi_0

%     InitialData_choice = 1; % Gaussian initial data
    InitialData_choice = 2; % Thomas Fermi initial data

%     X0 = 0; Y0 = 0; Z0 = 0;
%     gamma_x = 1;    gamma_y = 1;    gamma_z = 1;
    
    Phi_0 = InitialData_Var3d(Method, Geometry3D, Physics3D, InitialData_choice);
%     Phi_0 = InitialData_Var3d(Method, Geometry3D, Physics3D, InitialData_choice, X0, Y0, Z0, gamma_x, gamma_y, gamma_z);

    % introduce gaussian phase
    phase = GaussianInitialData3d(Geometry3D, Physics3D, 1, 1, 1, 0, 0, 0);
    
    for j = 1:Ncomponents
        Phi_0{j} = Phi_0{j} .* exp(1i*((j-1)*(2*pi/3) + phase));
    end

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
    globaluserdef_names{1} = 'Magnetization';
    globaluserdef_names{2} = 'Population psi+';
    globaluserdef_names{3} = 'Population psi0';
    globaluserdef_names{4} = 'Population psi-';
    
    Outputs = OutputsINI_Var3d(Method, Evo_outputs, Save_solution, [], [], ...
        globaluserdef_outputs,globaluserdef_names);

    %% Printing preliminary outputs

    Printing = 1;
    Evo = 50;
    Draw = 0;
    Print = Print_Var3d(Printing, Evo, Draw);

    %% RUN THE SIMULATION to find the ground state

    info.add_simulation_info(Geometry3D);
    [Phi_1, Outputs] = GPELab3d(Phi_0, Method, Geometry3D, Physics3D, Outputs, [], Print);
    
    save(info.get_workspace_path('phi_ini'), 'Phi_1')
    
    %% Save the workspace & simulation info

    % save information about final iteration in info file
    info.add_result_info(Method, Outputs);
    
    % saving groundstate workspace in v7.3 MAT file
    save(info.get_workspace_path('groundstate_v7.3'), '-v7.3');
    
    % saving groundstate workspace in 'normal' file
    save(info.get_workspace_path('groundstate'))

    %% Draw user-defined functions populations & magnetization
    
    close all;
    pause(2) % pauses the program for 2 seconds
    
    its = Outputs.Iterations;
    
    % Plot magnetization
    plot_magnetization(its, Outputs.User_defined_global{1}, info, Outputs.Evo_outputs)
    % Plot population fractions
    plot_populationfractions(its, Outputs.User_defined_global(2:4), info, Outputs.Evo_outputs)
    % Plot population distribution on x-axis
    plot_populationdistribution(Geometry3D, Phi_1, info)
    % Plot magnetization distribution on x-axis
    plot_magnetizationdistribution(Geometry3D, Phi_1, info)
    
    %% Time plots
    % magnetization distribution on x-axis
    timeslider_magnetizationdistribution(Geometry3D, Outputs.Solution, info)
    % population distribution on x-axis
    timeslider_populationdistribution(Geometry3D, Outputs.Solution, info)
    
    %% Draw & save solution

    close all;
    pause(1) % pauses the program for 1 second

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
            
    Draw_solution3d(Phi_0, Method, Geometry3D, Figure_Var3d());
    for i = 1 : Ncomponents
        info.save_figure(2*i-1, 'initialdata', phisqname{i},[],[]);
        info.save_figure(2*i, 'initialdata', anglename{i},[],[]);
        
        if ((Ncomponents > 1) && (i == Ncomponents))
            info.save_figure(2*i+1, 'initialdata', phisqname{i+1},[],[]);
        end
    end
    
    Draw_solution3d(Phi_1, Method, Geometry3D, Figure_Var3d());
    for i = 1 : Ncomponents
        info.save_figure(2*i-1, 'groundstate', phisqname{i},[],[]);
        info.save_figure(2*i, 'groundstate', anglename{i},[],[]);
        
        if ((Ncomponents > 1) && (i == Ncomponents))
            info.save_figure(2*i+1, 'groundstate', phisqname{i+1},[],[]);
        end
    end

    simulation_finished = 'Ground state simulation finished';

    save(info.get_workspace_path('fittingdata'),'simulation_finished', '-append')
    
    info.finish_info();
    
end