 %{
23Na/87Rb atom in the F=1 manifold
Spinor BEC in 3D
No rotations, no magnetic field, optical dipole + harmonic osc trap
%}

function [] = spinor_GPE3D_dynamics(info)
    
    close all;
    Phi_in = info.params.Phi_input;
   
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
    
%     %% Run a ground state simulation if Phi_in is not yet given.
%     
%     if ~exist('Phi_in', 'var')
%         groundst_simulation = 'spinor_GPE3D_ground';
%         ground_info = multi_runner(groundst_simulation, info.params.boxlimits, info.params.Ngridpts);
%         ws_ground = ground_info.get_workspace_path('groundstate'); % ground state workspace path/name
%         ws_ground = load(ws_ground); % ground state workspace struct of variables
%         names = fieldnames(ws_ground); % variable names cell in ground state ws
%         for i = 1:length(names)
%             if strcmp(names{i}, 'Phi_1')
%                 Phi_in = ws_ground.(names{21});
%             end
%         end
%         if ~exist('Phi_in', 'var')
%             error('Something went wrong with the ground state simulation, please check.')
%         end
%     end
        
    %{
    DYNAMICAL SIMULATION
    %}

    %% Simulation methods

    Computation = 'Dynamic';
    Ncomponents = 3;
    Type = 'Splitting'; % defaults to Strang splitting
%     Type = 'Relaxation';
    dx = (2*xlim / (Nx-1));
%     Deltat = info.params.dt;
    Deltat = min(0.25,0.1*(dx)^3);
    LimitingIter = 1000; % A limitation to the iterations bc time
    Stop_time = floor(min(100, (LimitingIter*Deltat)));
    Stop_crit = [];
    Max_iter = 1000;
%         Stop_crit = {'MaxNorm', 1e-12};
%         Precond_type = 'FLaplace'; % defaults to 'FLaplace'
%         Precond_type = []; % defaults to 'FLaplace'
%     Precond_type = 'FThomasFermi'; % defaults to 'FLaplace'
%     Output = 1;

    Method = Method_Var3d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit, Max_iter);
%     Method = Method_Var3d(Computation, Ncomponents, Type, Deltat, Stop_time, [], [], Precond_type, Output, Splitting);

    % Saving workspace with relevant data for fitting
    Method_dynamical = Method;
    Method.M = info.params.M;
    Method.projection = false;
    
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
    
    % Potential function
%     Udp0 = info.params.dipoleTrap0;
%     Wx = info.params.xOmega;
    gx = info.params.xOmega / info.params.trapfreq;
    gx = 10;
    gy = 1; gz = 1; gx = 1;
    
%     Wx = sqrt(4*Udp0 / (info.params.atom_mass*info.params.dipoleWaist_x^2));
%     Wy = getphysconst('hbar') / (info.params.atom_mass * getsimconst('axicon_radius')^2);
%     Wz = Wy;
%     Wx = min(Wx/100, Wy); % limiting this 2d behaviour
%     Wmin = min(min(Wx, Wy), Wz);
%     gx = Wx/Wmin; gy = Wy/Wmin; gz = Wz/Wmin; % scaled parameters gamma => 1
    
%     dipoletrap = dipoleTrap(info.params, X, Y, Z);
%     quadratictrap = quadratic_potential3d(gx, gy, gz, X, Y, Z);
    
%     Bz = 10^(-4); % Magnetic field in T (10^4 G = 1 T)
%     Bz = 10^(-3) * 10^(-4); % Magnetic field in T (10^4 G = 1 T)
%     Bz = 10^(-10); % 1 uG = 0.1nT = small field
%     Bz = 10^(-8); % 100 uG = 10 nT = moderate field
%     Bz = 10^(-4); % 1 G = 100 uT = big field
%     Bz = 100 * 10^(-4); % 100 G = 10 mT = even bigger field
    Bz = info.params.Bz;
%     Bz = 0;
%     Bz = 0; % no magnetic field
    
%     potential_with_Bfield = @(X,Y,Z) addingPotentials(info.params, ...
%         dipole_plus_quadratictrap(info.params, gx,gy,gz, X,Y,Z), ...
%         magneticFieldPotential(info.params, Bz, Wmin));
    % Optical dipole trap + harmonic oscillator (small-x approx) + axicon
    % in yz plane approximated as a harmonic oscillator with axicon radius
%     Nc = info.params.nComponents;
%     i = num2cell(eye(3));
%     fz = eye(Nc) .* [1, 0, -1]; fz2 = (fz).^2;
%     Fz = num2cell(fz); Fz2 = num2cell(fz2);
    
    % temp fix
%     gx = 1; gy = 1; gz = 1;
%     Wmin = 1 * 2*pi; % Hz
    p = info.params.p;
    q = info.params.q;
    
    % set q < 1/2
%     q = 0.25;
    % for a ferromagnetic states mix, need sqrt(2*q)<|p|<1 & 0<q<1/2
%     p = 0.5*(1+sqrt(2*q));
    
%     p = p*info.params.betas;
%     q = q*info.params.betas;
%     potential = cell(3,3);
%     for n = 1:Nc
%         for m = 1:Nc
%             if n == m
%                 potential{n,m} = @(X,Y,Z) ( i{n,m} * quadratic_potential3d(gx,gy,gz, X,Y,Z) ...
%                     + Fz2{n,m} * q);
% %                     + dipoleTrap(info.params, X,Y,Z) - Fz{n,m} * p + Fz2{n,m} * q);
%             end
%         end
%     end
    % making a simple potential
    potential = @(X,Y,Z) ( quadratic_potential3d(gx,gy,gz, X,Y,Z) );
                
%     pot = @(X,Y,Z) (eye(3) .* quadratic_potential3d(gx,gy,gz, X,Y,Z) + eye(3) .* dipoleTrap(info.params, X,Y,Z) ...
%         + (fx) * p + (fx).^2 * q );
%     pot = addingPotentialsAndBfield(info.params, ...
%         @(X,Y,Z) (@(X,Y,Z) quadratic_potential3d(gx,gy,gz, X,Y,Z) + @(X,Y,Z) dipoleTrap(info.params, X,Y,Z)), ...
%         @(l) magneticFieldPotential(info.params, Bz, Wmin, l));
    
%     Physics3D = Potential_Var3d(Method, Physics3D, potential_with_Bfield, ones(Nc,Nc));
%     Physics3D = Potential_Var3d(Method, Physics3D, potential, ones(Nc,Nc));
    Physics3D = Potential_Var3d(Method, Physics3D, potential);

    % Nonlinearity
%     Physics3D = Nonlinearity_Var3d(Method, Physics3D, Coupled_Cubic3d_spin1(Betan,Betas)); % cubic nonlinearity with off-diagonal coupling
    Physics3D = Nonlinearity_Var3d(Method, Physics3D, Coupled_Cubic3d_spin1(Betan,Betas,p,q), [], ...
        Coupled_CubicEnergy3d_spin1(Betan,Betas,p,q)); % cubic nonlinearity with off-diagonal coupling
    
    %% version mgmt
    curdir = strsplit(pwd, '/');
    curdir = strsplit(curdir{end}, '\');
    curdir = curdir{end};

    %% Printing interaction strength

    aho = sqrt(getphysconst('hbar') / (info.params.atom_mass * info.params.trapfreq));
    
    info.add_info_separator();
    info.add_custom_info('Folder: \t %s \n', curdir); % print current folder
    info.add_info_separator();
    info.add_custom_info('atom \t=\t %s \n', info.params.atom);
    info.add_custom_info('Beta_n \t=\t %f \n', info.params.betan); % print interaction parameter Beta_n
    info.add_custom_info('Beta_s \t=\t %f \n', info.params.betas); % print interaction parameter Beta_s
    info.add_custom_info('Wmin \t=\t %.1g x 2pi Hz \n', info.params.trapfreq / (2*pi)); % (minimum) trap frequency
    info.add_custom_info('gammas \t=\t [%.4f,%.4f,%.4f] \n', gx,gy,gz); % print gammas
    info.add_custom_info('Magnetic field parameters:\n'); % 
    info.add_custom_info('\tBz \t=\t %.3g Gauss\n', Bz*10^4); % magnetic field in Gauss
    info.add_custom_info('\tLinear and quadratic Zeeman energy,\n');
    info.add_custom_info('\t[p,q] \t=\t %.3g, %.3g (in h.o. energy units)\n', p,q); % Zeeman energy
    info.add_custom_info('\t \t=\t %.3g, %.3g (in units beta_s)\n', p/info.params.betas,q/info.params.betas); % Zeeman energy
    info.add_custom_info('dt \t=\t %f \n', info.params.dt);
    info.add_custom_info('dx \t=\t %f \n', dx);
    info.add_info_separator();
    info.add_custom_info('a_ho \t=\t %f um\n', aho*10^(6)); % harmonic oscillator length in um
    info.add_custom_info('a_ho(xyz) \t=\t [%.2g,%.2g,%.2g] aho, or:\n', ...
        1/sqrt(gx),1/sqrt(gy),1/sqrt(gz)); % normalized harm osc length (x,y,z)
    info.add_custom_info('a_ho(xyz) \t=\t [%.2g,%.2g,%.2g] um, or:\n', ...
        aho*10^(6)/sqrt(gx), aho*10^(6)/sqrt(gy), aho*10^(6)/sqrt(gz)); % harm osc length (x,y,z)
    info.add_info_separator();
    
    %% Determining outputs

    % Must be equal to or smaller than Evo from Print
    Evolim = round((3*(Nx)^3*Stop_time / (Deltat*7e7))/5)*5;
    Evo_outputs = max(10, Evolim);
    if Max_iter < 101
        Evo_outputs = min(5,Evo_outputs);
    end
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
    globaluserdef_outputs{8} = @(Phi,X,Y,Z,FFTX,FFTY,FFTZ) ...
        Directional_Magnetization(Method, Geometry3D, Phi, 'F2', X, Y, Z, FFTX, FFTY, FFTZ);
    globaluserdef_outputs{9} = @(Phi,X,Y,Z,FFTX,FFTY,FFTZ) ...
        Directional_Magnetization(Method, Geometry3D, Phi, 'M2', X, Y, Z, FFTX, FFTY, FFTZ);
    globaluserdef_names{1} = 'Longitudinal magnetization';
    globaluserdef_names{2} = 'Population psi+';
    globaluserdef_names{3} = 'Population psi0';
    globaluserdef_names{4} = 'Population psi-';
    globaluserdef_names{5} = 'Magnetization Fx';
    globaluserdef_names{6} = 'Magnetization Fy';
    globaluserdef_names{7} = 'Magnetization Fz';
    globaluserdef_names{8} = 'Total magnetization F^2';
    globaluserdef_names{9} = 'Total magnetization M^2';
    
    Outputs = OutputsINI_Var3d(Method, Evo_outputs, Save_solution, [], [], ...
        globaluserdef_outputs,globaluserdef_names);
    
    %% Printing preliminary outputs
    
    Printing = 1;
    % Must be equal to or bigger than Evo_outputs
    Evo = max(25, Evo_outputs);
    if Max_iter < 60
        Evo = 10;
    end
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
    F = Outputs.User_defined_global(5:7);
    
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
    % phi distribution as a sliced 3d function
    timeslider_slicer(Geometry3D, Outputs.Solution, info)
    
    % adding extra parameters to workspace
    save(info.get_workspace_path('dynamics_v7.3'), 'F', 'its', '-append');
    
    close all;
    sprintf('%s', info.get_workspace_path('dynamics_v7.3'))
    
    %% Save PhiData structures for fitting w/o the whole workspace

    phi_input = PhiData(Phi_in, Geometry3D);
    phi_ground = PhiData(Phi_in, Geometry3D);
    phi_dyn = PhiData(Phi, Geometry3D);
    timeunits = 1/info.params.trapfreq;
    
    % Saving workspace with relevant data for fitting
    save(info.get_workspace_path('fittingdata'), 'gx', 'gy', 'gz', ...
        'phi_dyn', 'phi_ground', 'phi_input', 'info', 'Geometry3D', ...
        'Method', 'timeunits', '-append'); % to not overwrite Method_ground

    simulation_finished = 'yes!';
    save(info.get_workspace_path('fittingdata'),'simulation_finished', '-append')
    
    info.finish_info();
    
    %% end
end
