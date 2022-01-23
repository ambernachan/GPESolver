 %{
23Na/87Rb atom in the F=1 manifold
Spinor BEC in 1D
No rotations, no magnetic field, optical dipole + harmonic osc trap
%}

function [] = spinor_GPE1D_dynamics(info)
    
    close all;
    
    Phi_in = info.params.Phi_input;
   
    %% Setting variables
    
    % Setting simulation space
    xlim = info.params.boxlimits(1);
    Nx = info.params.Ngridpts;
    
    % Setting physical parameters
    if isprop(info.params, 'delta')
        Delta = info.params.delta;
    else
        Delta = 0.5;
        info.params.delta = Delta;
    end
    if isprop(info.params, 'gammas')
        gx = info.params.gammas(1);
    else
        gx = 1;
        info.params.gammas(1) = gx;
    end
    if isprop(info.params, 'dimensions')
        dimensions = info.params.dimensions;
    else
        dimensions = 1;
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
%         groundst_simulation = 'spinor_GPE1D_ground';
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
    if dx >= 1
        warning('Your simulation may fail because dx => 1.')
    end
%     Deltat = info.params.dt;
%     Deltat = min(min(min(0.0625,0.1*(dx)^3),0.00625),0.0000625);
    Deltat = info.params.dt; 
    if ~isempty(info.params.iterations)
        LimitingIter = info.params.iterations;
        Max_iter = info.params.iterations;
    else
        LimitingIter = 10000; % A limitation to the iterations bc time
        Max_iter = 10000;
    end
%     Stop_time = floor(min(1000, (LimitingIter*real(Deltat))));
    Stop_time = max(1,floor(LimitingIter*real(Deltat)));
    Stop_crit = [];
    
%     Stop_time = floor(min(1000, round(Max_iter*Deltat/5)*5));
    
%         Stop_crit = {'MaxNorm', 1e-12};
%         Precond_type = 'FLaplace'; % defaults to 'FLaplace'
%         Precond_type = []; % defaults to 'FLaplace'
%     Precond_type = 'FThomasFermi'; % defaults to 'FLaplace'
%     Output = 1;

    Method = Method_Var1d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit, Max_iter);
%     Method = Method_Var1d(Computation, Ncomponents, Type, Deltat, Stop_time, [], [], Precond_type, Output, Splitting);

    % Saving workspace with relevant data for fitting
    Method_dynamical = Method;
    Method.M = info.params.M;
    Method.projection = info.params.projection;
    Method.q = info.params.q; % Required for stability of the simulation!
    
    %% Save dynamical method to fittingdata workspace
    if isfile(info.get_workspace_path('fittingdata'))
        save(info.get_workspace_path('fittingdata'), 'Method_dynamical', '-append');
    else
        save(info.get_workspace_path('fittingdata'), 'Method_dynamical');
    end
    clear Method_dynamical;

    %% Geometry1D
    
    xmin = -xlim;   xmax = xlim;
    
    Geometry1D = Geometry1D_Var1d(xmin, xmax, Nx);%% Geometry1D
        
    %% Physics1D

    % Delta is already defined
    Beta = info.params.beta; % multiplication factor for Beta_n and Beta_s
    Betan = info.params.betan;
    Betas = info.params.betas;
    Physics1D = Physics1D_Var1d(Method, Delta, Beta);
    Physics1D = Dispersion_Var1d(Method, Physics1D); % !!!
    
    % Potential function
    gx = info.params.xOmega / info.params.trapfreq;
    gx = 1;
    Bz = info.params.Bz;
    Bmin = info.params.Bmin;
    
    p = info.params.p;
    q = info.params.q;
    
    % making a simple potential
    potential = @(X) ( quadratic_potential1d(gx, X) );
                
%     Physics1D = Potential_Var1d(Method, Physics1D, potential_with_Bfield, ones(Nc,Nc));
%     Physics1D = Potential_Var1d(Method, Physics1D, potential, ones(Nc,Nc));
    Physics1D = Potential_Var1d(Method, Physics1D, potential);

    % Nonlinearity
%     Physics1D = Nonlinearity_Var1d(Method, Physics1D, Coupled_Cubic1d_spin1(Betan,Betas)); % cubic nonlinearity with off-diagonal coupling
    Physics1D = Nonlinearity_Var1d(Method, Physics1D, Coupled_Cubic1d_spin1(Betan,Betas,info.params), [], ...
        Coupled_CubicEnergy1d_spin1(Betan,Betas));%,info.params)); % cubic nonlinearity with off-diagonal coupling
    
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
    info.add_custom_info('gamma \t=\t [%.4f,%.4f,%.4f] \n', gx); % print gamma
    info.add_custom_info('Magnetic field parameters:\n'); % 
    info.add_custom_info('\tBz \t=\t %.3g Gauss\n', Bz*10^4); % magnetic field in Gauss
    info.add_custom_info('\tLinear and quadratic Zeeman energy,\n');
    if Bmin == 0
        info.add_custom_info('\t[p,q] \t=\t %.3g, %.3g (in h.o. energy units)\n', p,q); % Zeeman energy
        info.add_custom_info('\t \t=\t %.3g, %.3g (in units beta_s)\n', p/info.params.betas,q/info.params.betas); % Zeeman energy
    else
        [pmin,qmin] = getMagneticFieldPars(Bmin, info.params.trapfreq, info.params.Ehfs);
        info.add_custom_info('\t[p,q] \t=\t (%.3g;%.3g), (%.3g;%.3g) (in h.o. energy units)\n', pmin,p,qmin,q); % Zeeman energy
        info.add_custom_info('\t \t=\t (%.3g;%.3g), (%.3g;%.3g) (in units beta_s)\n', ...
            pmin/info.params.betas,p/info.params.betas,qmin/info.params.betas,q/info.params.betas); % Zeeman energy
    end
    sgn_imag_dt = sign(imag(Deltat)); str_sgn = '+'; 
    if (sgn_imag_dt == -1) str_sgn = '-'; end
    info.add_custom_info('dt \t=\t %.2g%s%.2gi \n', real(Deltat), str_sgn, abs(imag(Deltat)));
    info.add_custom_info('\t (dt = %.3g us real time & %.2g pct damping) \n', ...
        real(Deltat)*1e6/info.params.trapfreq, abs(imag(Deltat)) / abs(Deltat));
%     info.add_custom_info('dt \t=\t %f \n', info.params.dt);
    info.add_custom_info('dx \t=\t %f \n', dx);
    info.add_info_separator();
    info.add_custom_info('a_ho \t=\t %.2g um\n', aho*10^(6)); % harmonic oscillator length in um
    info.add_custom_info('a_ho(x) \t=\t %.2g aho, or:\n', ...
        1/sqrt(gx)); % normalized harm osc length (x)
    info.add_custom_info('a_ho(x) \t=\t %.2g um, or:\n', ...
        aho*10^(6)/sqrt(gx)); % harm osc length (x)
    info.add_info_separator();
    
    %% Determining outputs

    % Must be equal to or smaller than Evo from Print
    Evolim = round(((3*(Nx)^3*Stop_time / (real(Deltat)*7e7))^(1/3))/5)*5;
    Evo_outputs = max(10, Evolim);
    if Max_iter < 101
        Evo_outputs = min(5,Evo_outputs);
    end
    Save_solution = 1;
    
    globaluserdef_outputs{1} = @(Phi,X,FFTX) ...
        Magnetization(Method, Geometry1D, Phi, X,[],[], FFTX,[],[]);
    globaluserdef_outputs{2} = @(Phi,X,FFTX) ...
        Population(Method, Geometry1D, Phi, X,[],[], FFTX,[],[], 1);
    globaluserdef_outputs{3} = @(Phi,X,FFTX) ...
        Population(Method, Geometry1D, Phi, X,[],[], FFTX,[],[], 0);
    globaluserdef_outputs{4} = @(Phi,X,FFTX) ...
        Population(Method, Geometry1D, Phi, X,[],[], FFTX,[],[], -1);
    globaluserdef_outputs{5} = @(Phi,X,FFTX) ...
        Directional_Magnetization(Method, Geometry1D, Phi, 'x', X,[],[], FFTX,[],[]);
    globaluserdef_outputs{6} = @(Phi,X,FFTX) ...
        Directional_Magnetization(Method, Geometry1D, Phi, 'y', X,[],[], FFTX,[],[]);
    globaluserdef_outputs{7} = @(Phi,X,FFTX) ...
        Directional_Magnetization(Method, Geometry1D, Phi, 'z', X,[],[], FFTX,[],[]);
    globaluserdef_outputs{8} = @(Phi,X,FFTX) ...
        Directional_Magnetization(Method, Geometry1D, Phi, 'M2', X,[],[], FFTX,[],[]);
    globaluserdef_outputs{9} = @(Phi,X,FFTX) ...
        Directional_Magnetization(Method, Geometry1D, Phi, 'F2', X,[],[], FFTX,[],[]);
    globaluserdef_names{1} = 'Longitudinal magnetization';
    globaluserdef_names{2} = 'Population psi+';
    globaluserdef_names{3} = 'Population psi0';
    globaluserdef_names{4} = 'Population psi-';
    globaluserdef_names{5} = 'Magnetization Fx';
    globaluserdef_names{6} = 'Magnetization Fy';
    globaluserdef_names{7} = 'Magnetization Fz';
    globaluserdef_names{8} = 'Total magnetization M^2';
    globaluserdef_names{9} = 'Total magnetization F^2';
    
    Outputs = OutputsINI_Var1d(Method, Evo_outputs, Save_solution, [], [], ...
        globaluserdef_outputs,globaluserdef_names);
    
    %% Printing preliminary outputs
    
    Printing = 0;
    % Must be equal to or bigger than Evo_outputs
    Evo = max(25, Evo_outputs);
    if Max_iter < 60
        Evo = 10;
    end
    Draw = 0;
    Print = Print_Var1d(Printing, Evo, Draw);

    %% RUN THE SIMULATION to find the ground state

    info.add_simulation_info(Geometry1D);
%         [Phi, Outputs] = GPELab1d(Phi_i, Method, Geometry1D, Physics1D, Outputs, [], Print);
    [Phi, Outputs] = GPELab1d(Phi_in, Method, Geometry1D, Physics1D, Outputs, [], Print);
    
    save(info.get_workspace_path('phi_ini'), 'Phi_in', 'Phi', 'Outputs', 'info', '-v7.3')
    
    %% Save the workspace & simulation info

    % save information about final simulation iteration in info file
    info.add_result_info(Method, Outputs);

    % saving dynamics workspace in v7.3 MAT file
    save(info.get_workspace_path('dynamics_v7.3'), '-v7.3');
    
    % saving dynamics workspace in 'normal' file
    save(info.get_workspace_path('dynamics'))
    
    %% create file that shows type of atom and simulation stats
    fname = createTextFileName(info, Geometry1D, Method, Outputs.Iterations*Outputs.Evo_outputs);
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
    plot_populationdistribution(Geometry1D, Phi, info, 'x')
    
    hold on;
    BZ = @(z) info.params.Bmin + (info.params.Bz-info.params.Bmin)*(1+z/xlim)/2;
    zz = -xlim:dx:xlim;
    yyaxis right
    ylabel('Bz (G)')
    plot(zz,BZ(zz)*10^4)
    ax2 = gca;
    ax2.Children.DisplayName = 'Magnetic field';
    
    % Plot magnetization distribution on x-axis
    plot_magnetizationdistribution(Geometry1D, Phi, info)
    % Plot transverse & longitudinal magnetization
    plot_magnetizations(its, F, info, Outputs.Evo_outputs, Method)
    
    %% Time plots
    % magnetization distribution on x-axis
    timeslider_magnetizationdistribution(Geometry1D, Outputs.Solution, info)
    % population distribution on x-axis
    timeslider_populationdistribution(Geometry1D, Outputs.Solution, info)
    % phase distribution as a sliced 1d function
    timeslider_phase(Geometry1D, Outputs.Solution, info)
    % phi distribution as a sliced 1d function
    timeslider_slicer(Geometry1D, Outputs.Solution, info)
    
    % Gaussian + Thomas-Fermi comparison
    compareGTF(Geometry1D, Phi, info);
    
    % adding extra parameters to workspace
    save(info.get_workspace_path('dynamics_v7.3'), 'F', 'its', '-append');
    
    close all;
    sprintf('%s', info.get_workspace_path('dynamics_v7.3'))
    
    %% Save PhiData structures for fitting w/o the whole workspace

    phi_input = PhiData(Phi_in, Geometry1D);
    phi_ground = PhiData(Phi_in, Geometry1D);
    phi_dyn = PhiData(Phi, Geometry1D);
    timeunits = 1/info.params.trapfreq;
    
    % Saving workspace with relevant data for fitting
    save(info.get_workspace_path('fittingdata'), 'gx', ...
        'phi_dyn', 'phi_ground', 'phi_input', 'info', 'Geometry1D', ...
        'Method', 'timeunits', '-append'); % to not overwrite Method_ground

    simulation_finished = 'yes!';
    save(info.get_workspace_path('fittingdata'),'simulation_finished', '-append')
    
    info.finish_info();
    
    %% end
end
