%{
23Na in the F=1 manifold
Spinor BEC in 3D
No rotations, no magnetic field, quadratic optical trap
%}

function [] = spinor_GPE3D_ground(info)
    
    close all;
   
    %% Setting variables
    % Determine interaction strength compared to kinetic energy
    if isfield(info.params, 'a0')
        a0 = info.params.a0;
    else
        a0 = getsimconst('a0_Na');
        info.params.a0 = a0;
    end
    if isfield(info.params, 'a2')
        a2 = info.params.a2;
    else
        a2 = getsimconst('a2_Na');
        info.params.a2 = a2;
    end
    N = getsimconst('N'); % number of particles
    hbar = getphysconst('hbar'); % in kg m^2 / s
    trapfreq = getsimconst('trap_freq'); % Trap strength in Hz (symmetric for now)
    atom_mass = getsimconst('mass_Na'); % Atom mass in kg
    spin_pair = getsimconst('spin_pair'); % hyperfine spin manifold (=1)
    info.params.spin_pair = spin_pair;
    
    aho = sqrt(hbar / (atom_mass * trapfreq));
    chin = N*(2*a2+a0)/(3*aho);
    chis = N*(a2-a0)/(3*aho);
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
    
    %%
    %{
    GROUND STATE SIMULATION TO GET STARTING FUNCTION PHI_0 FOR DYNAMIC SIMULATION
    %}

    %% Simulation methods

    Computation = 'Ground';
    Ncomponents = 3;
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
    %Physics3D = Potential_Var3d(Method, Physics3D); % std quadratic potential
    %Physics3D = Potential_Var3d(Method, Physics3D, @(X,Y,Z) quadratic_potential3d(1, 5, 1, X, Y, Z));
    Physics3D = Potential_Var3d(Method, Physics3D, @(X,Y,Z) quadratic_potential3d(gx, gy, gz, X, Y, Z));
    
    Physics3D = Nonlinearity_Var3d(Method, Physics3D, Coupled_Cubic3d_spin1(Betan,Betas), [], ...
        Coupled_CubicEnergy3d_spin1(Betan,Betas)); % cubic nonlinearity with off-diagonal coupling
    
    %Physics3D = Gradientx_Var3d(Method, Physics3D, @(x,y) -1i*Omega*y);
    %Physics3D = Gradienty_Var3d(Method, Physics3D, @(x,y) 1i*Omega*x);

    %% Defining a starting function Phi_0

    InitialData_choice = 1; % Gaussian initial data

%     % find w3d
%     exp = ExpPhis(0, [gx,gy,gz], [0,0,0], Geometry3D);
%     W = exp.getW();
%     w = W.d3;
%     clear exp W
% 
%     sigma_w_delta = struct();
%     sigma_w_delta.x = w * (Delta/2)^(0.25) / sqrt(gx); sigma_w_delta.x = sigma_w_delta.x(sigma_w_delta.x>0);
%     sigma_w_delta.y = w * (Delta/2)^(0.25) / sqrt(gy); sigma_w_delta.y = sigma_w_delta.y(sigma_w_delta.y>0);
%     sigma_w_delta.z = w * (Delta/2)^(0.25) / sqrt(gz); sigma_w_delta.z = sigma_w_delta.z(sigma_w_delta.z>0);

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
%     info.add_custom_info('S \t=\t %f \n', 0); % print interaction strength
    info.add_custom_info('Beta_n \t=\t %f \n', Betan); % print interaction parameter Beta_n
    info.add_custom_info('Beta_s \t=\t %f \n', Betas); % print interaction parameter Beta_s
%     info.add_custom_info('w \t=\t %f \n', w); % print Gaussian parameter w
%     info.add_custom_info('sigmas \t=\t [%.4f,%.4f,%.4f] \n', ...
%         sigma_w_delta.x,sigma_w_delta.y,sigma_w_delta.z); % print Gaussian parameter sigma_xyz with Delta
    info.add_custom_info('Delta \t=\t %f \n', Delta); % print kinetic energy parameter Delta
    info.add_custom_info('gammas \t=\t [%.4f,%.4f,%.4f] \n', gx,gy,gz); % print gammas
    info.add_info_separator();

    %% Determining outputs
    
    Save_solution = 1;
    Outputs = OutputsINI_Var3d(Method, Save_solution);

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
    
end