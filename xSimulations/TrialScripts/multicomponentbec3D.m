%{
Multi-component system without rotations in 3D
%}

function [] = multicomponentbec3D(info)

    close all;

    %% Setting variables
    % Determine interaction strength compared to kinetic energy
    S = info.params.S;
    % Setting simulation space
    xlim = info.params.boxlimits(1); ylim = info.params.boxlimits(2); zlim = info.params.boxlimits(3);
    %Nx = Ngridpts(1); Ny = Ngridpts(2); Nz = Ngridpts(3);
    Nx = info.params.Ngridpts; Ny = info.params.Ngridpts; Nz = info.params.Ngridpts;
    % Setting physical parameters
    if isfield(info.params, 'delta')
        Delta = info.params.delta;
    else
        Delta = 0.5;
    end
    gx = info.params.gammas(1); gy = info.params.gammas(2); gz = info.params.gammas(3);
    
    %%
    %{
    GROUND STATE SIMULATION
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
    
    xmin = -xlim; xmax = xlim;
    ymin = -ylim; ymax = ylim;
    zmin = -zlim; zmax = zlim;
    
    Geometry3D = Geometry3D_Var3d(xmin, xmax, ymin, ymax, zmin, zmax, Nx, Ny, Nz);

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
    

end
