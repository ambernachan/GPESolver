%{
23Na in the F=1 manifold
Spinor BEC in 1D
No rotations, no magnetic field, quadratic optical trap
One can change the type of atom by changing :
    - the atom mass
    - the a0 and a2 parameters
%}

function [] = spinor_GPE1D_ground(info)
    
    close all;
   
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
    
    %%
    %{
    GROUND STATE SIMULATION TO GET STARTING FUNCTION PHI_0 FOR DYNAMIC SIMULATION
    %}

    %% Simulation methods

    Computation = 'Ground';
    Ncomponents = 3;
    Type = 'BESP';
    dx = (2*xlim / (Nx-1));
    if dx >= 1
        warning('You simulation may fail because dx => 1.')
    end
%     Deltat = 0.1*dx^3;
    Deltat = 0.0625;
    Stop_crit = {'MaxNorm', 1e-50};
    Max_iter = 137500;
    Stop_time = floor(min(1000, round(Max_iter*Deltat/5)*5));
%     LimitingIter = 1000; % A limitation to the iterations bc time
%     Stop_time = floor(min(100, (LimitingIter*Deltat)));
%     Stop_crit = [];
%     Max_iter = 1000;

    Method = Method_Var1d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit, Max_iter);
    Method.M = info.params.M;
    Method.q = info.params.q; % Required for stability of the simulation!
    Method.projection = info.params.projection;

    %% Geometry1D
    
    xmin = -xlim;   xmax = xlim;
    
    Geometry1D = Geometry1D_Var1d(xmin, xmax, Nx);

    %% Take into account possible higher # of chis/S
    
%     Save ground state method
    Method_ground = Method;
    save(info.get_workspace_path('fittingdata'),'Method_ground')
    clear Method_ground;
    
    %% Physics1D

    % Delta is already defined
    Beta = info.params.beta; % multiplication factor for Beta_n and Beta_s
    Betan = info.params.betan;
    Betas = info.params.betas;
    Physics1D = Physics1D_Var1d(Method, Delta, Beta);
    Physics1D = Dispersion_Var1d(Method, Physics1D); % !!!
    
    gx = info.params.xOmega / info.params.trapfreq;
    gx = 1;
    Bz = info.params.Bz;
    Bmin = info.params.Bmin;
    
    p = info.params.p;
    q = info.params.q;
    
    % making a simple potential
    potential = @(X) ( quadratic_potential1d(gx, X) );
    
    Physics1D = Potential_Var1d(Method, Physics1D, potential);

    % Nonlinearity
%     Physics1D = Nonlinearity_Var1d(Method, Physics1D); % std cubic nonlinearity
%       Physics1D = Nonlinearity_Var1d(Method, Physics1D, Coupled_Cubic1d_spin1(Betan,Betas)); % cubic nonlinearity with off-diagonal coupling
    Physics1D = Nonlinearity_Var1d(Method, Physics1D, Coupled_Cubic1d_spin1(Betan,Betas, info.params), [], ...
        Coupled_CubicEnergy1d_spin1(Betan,Betas, info.params)); % cubic nonlinearity with off-diagonal coupling

    %% Defining a starting function Phi_0

%     InitialData_choice = 1; % Gaussian initial data
%     InitialData_choice = 2; % Thomas Fermi initial data
%     Phi_0 = InitialData_Var1d(Method, Geometry1D, Physics1D, InitialData_choice);

    Phi_0 = info.params.Phi_input;

    %% version mgmt
    curdir = strsplit(pwd, '/');
    curdir = strsplit(curdir{end}, '\');
    curdir = curdir{end};

    %% Printing interaction strength
    
    aho = sqrt(getphysconst('hbar') / (info.params.atom_mass * info.params.trapfreq));
    avg_density = theDensity(Phi_0, Geometry1D);
    Us = info.params.betas * avg_density;
    
    info.add_info_separator();
    info.add_custom_info('Folder: \t %s \n', curdir); % print current folder
    info.add_info_separator();
    info.add_custom_info('atom \t=\t %s \n', info.params.atom);
    info.add_custom_info('Beta_n \t=\t %f \n', info.params.betan); % print interaction parameter Beta_n
    info.add_custom_info('Beta_s \t=\t %f \n', info.params.betas); % print interaction parameter Beta_s
    info.add_custom_info('Interaction energies,\n\t[Bn,Bs]*||phi||^2 \t=\t %.3g, %.3g\n', ...
        info.params.betan*avg_density,info.params.betas*avg_density); % interaction energy
    info.add_custom_info('Wmin \t=\t %.1g x 2pi Hz \n', info.params.trapfreq / (2*pi)); % (minimum) trap frequency
    info.add_custom_info('gamma \t=\t %.4f \n', gx); % print gamma
    info.add_custom_info('Magnetic field parameters:\n'); % 
    if Bmin == 0
        info.add_custom_info('\tBz \t=\t %.3g Gauss\n', Bz*10^4); % magnetic field in Gauss
        info.add_custom_info('\tLinear and quadratic Zeeman energy,\n');
        info.add_custom_info('\t[p,q] \t=\t %.3g, %.3g (in h.o. energy units)\n', p,q); % Zeeman energy
        info.add_custom_info('\t \t=\t %.3g, %.3g (in units beta_s*||phi||^2)\n', p/Us,q/Us); % Zeeman energy
    else
        info.add_custom_info('\tBz \t=\t [%.3g, %.3g] Gauss\n', Bmin*10^4, Bz*10^4); % magnetic field in Gauss
        info.add_custom_info('\tLinear and quadratic Zeeman energy,\n');
        [pmin,qmin] = getMagneticFieldPars(info.params.Bmin, info.params.trapfreq, info.params.Ehfs);
        info.add_custom_info('\t[p,q] \t=\t (%.3g;%.3g), (%.3g;%.3g) (in h.o. energy units)\n', pmin,p,qmin,q); % Zeeman energy
        info.add_custom_info('\t \t=\t (%.3g;%.3g), (%.3g;%.3g) (in units beta_s*||phi||^2)\n', ...
            pmin/Us,p/Us, qmin/Us,q/Us); % Zeeman energy
    end
    if Method.projection proj = 'yes'; else proj = 'no'; end
    info.add_custom_info('projections used? [%s] \t M = %.1g \n', proj, Method.M); % tells user whether projection constants are implemented
    info.add_custom_info('dt \t=\t %f \n', Deltat);
    info.add_custom_info('dx \t=\t %f \n', dx);
    info.add_info_separator();
    info.add_custom_info('a_ho \t=\t %.2g um\n', aho*10^(6)); % harmonic oscillator length in um
    info.add_custom_info('a_ho(x) \t=\t %.2g aho, or:\n', ...
        1/sqrt(gx)); % normalized harm osc length (x,y,z)
    info.add_custom_info('a_ho(x) \t=\t %.2g um, or:\n', ...
        aho*10^(6)/sqrt(gx)); % harm osc length (x,y,z)
    info.add_info_separator();

    %% Determining outputs
    
    % Must be equal to or smaller than Evo from Print
    Evolim = round(((3*(Nx)^3*Stop_time / (Deltat*7e7))^(1/3))/5)*5;
    Evo_outputs = max(10, Evolim);
    if Max_iter < 101
        Evo_outputs = min(5,Evo_outputs);
    end
%     Evo_outputs = 1;
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

    Printing = 1;
    % Must be equal to or bigger than Evo_outputs
    Evo = max(25, Evo_outputs);
    if Max_iter < 60
        Evo = 10;
    end
%     Evo = 1;
    Draw = 0;
    Print = Print_Var1d(Printing, Evo, Draw);

    %% RUN THE SIMULATION to find the ground state

    info.add_simulation_info(Geometry1D);
    [Phi_1, Outputs] = GPELab1d(Phi_0, Method, Geometry1D, Physics1D, Outputs, [], Print);
    
    save(info.get_workspace_path('phi_ini'), 'Phi_1')
    
    avg_density = theDensity(Phi_1, Geometry1D);
    Us = info.params.betas * avg_density;
    info.add_info_separator();
    info.add_custom_info('\tInteraction energies,\n\t[Bn,Bs]*||phi||^2 \t=\t %.3g, %.3g\n', ...
        info.params.betan*avg_density,info.params.betas*avg_density); % interaction energy
    if Bmin == 0
        info.add_custom_info('Linear and quadratic Zeeman energy,\n');
        info.add_custom_info('\t[p,q] \t=\t %.3g, %.3g (in h.o. energy units)\n', p,q); % Zeeman energy
        info.add_custom_info('\t \t=\t %.3g, %.3g (in units beta_s*||phi||^2)\n', p/Us,q/Us); % Zeeman energy
    else
        info.add_custom_info('Linear and quadratic Zeeman energy,\n');
        [pmin,qmin] = getMagneticFieldPars(info.params.Bmin, info.params.trapfreq, info.params.Ehfs);
        info.add_custom_info('\t[p,q] \t=\t (%.3g;%.3g), (%.3g;%.3g) (in h.o. energy units)\n', pmin,p,qmin,q); % Zeeman energy
        info.add_custom_info('\t \t=\t (%.3g;%.3g), (%.3g;%.3g) (in units beta_s*||phi||^2)\n', ...
            pmin/Us,p/Us, qmin/Us,q/Us); % Zeeman energy
    end
    info.add_info_separator();
    
    %% Save the workspace & simulation info

    % save information about final iteration in info file
    info.add_result_info(Method, Outputs);
    
    % saving groundstate workspace in v7.3 MAT file
    save(info.get_workspace_path('groundstate_v7.3'), '-v7.3');
    
    % saving groundstate workspace in 'normal' file
    save(info.get_workspace_path('groundstate'))

    %% create file that shows type of atom and simulation stats
    fname = createTextFileName(info, Geometry1D, Method, Outputs.Iterations*Outputs.Evo_outputs);
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
    plot_populationfractions(its, Outputs.User_defined_global(2:4), info, Outputs.Evo_outputs)
    % Plot population distribution on x-axis
    plot_populationdistribution(Geometry1D, Phi_1, info, 'x')
    
    hold on;
    BZ = @(z) info.params.Bmin + (info.params.Bz-info.params.Bmin)*(1+z/xlim)/2;
    zz = -xlim:dx:xlim;
    yyaxis right
    ylabel('Bz (G)')
    plot(zz,BZ(zz)*10^4)
    ax2 = gca;
    ax2.Children.DisplayName = 'Magnetic field';
        % Save figure
        savename = 'Population distribution w Bfield';
        info.save_figure(1, savename, '')
        info.save_figure(1, savename, '', info.fulldir, '.png')
        hold off
    
    % Plot magnetization distribution on z-axis
    plot_magnetizationdistribution(Geometry1D, Phi_1, info, 'x')
    % Plot transverse & longitudinal magnetization
    plot_magnetizations(its, F, info, Outputs.Evo_outputs, Method)
    
    %% Time plots
    % magnetization distribution on z-axis
    timeslider_magnetizationdistribution(Geometry1D, Outputs.Solution, info, 'x')
    % population distribution on x- and z-axis
    timeslider_populationdistribution(Geometry1D, Outputs.Solution, info, 'x')
    timeslider_populationdistribution(Geometry1D, Outputs.Solution, info, 'x')
    % phase distribution as a sliced 1d function
    timeslider_phase(Geometry1D, Outputs.Solution, info)
    % phi distribution as a sliced 1d function
    timeslider_slicer(Geometry1D, Outputs.Solution, info)
    
    % Gaussian + Thomas-Fermi comparison
    compareGTF(Geometry1D, Phi_1, info);
    
    % saving groundstate workspace in v7.3 MAT file
    save(info.get_workspace_path('groundstate_v7.3'), 'F', 'its', '-append');
    
    %% Draw & save solution

    close all;
    sprintf('%s', info.get_workspace_path('groundstate_v7.3'))

    simulation_finished = 'Ground state simulation finished';
    save(info.get_workspace_path('fittingdata'),'simulation_finished', '-append')
    
    info.finish_info();
    
end