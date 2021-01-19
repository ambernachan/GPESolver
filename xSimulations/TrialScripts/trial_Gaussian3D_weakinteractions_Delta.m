%{
...
%}

function [] = trial_Gaussian3D_weakinteractions_Delta(chi,xlimit,ylimit,zlimit,...
    xparticles,yparticles,zparticles,delta,gammas,branch)

%     clear;
%     chi = 0.01; xlimit = 3; xparticles = 2^5+1; ylimit = xlimit; yparticles = 2^5+1;
    close all;
    clearvars -except chi xlimit ylimit zlimit xparticles yparticles zparticles delta gammas branch
    %clear;

    %% Determine interaction strength compared to kinetic energy

    S = chi;

    %% Saving info files

    dimensions = 3;
    info = Info(name_from_filename(mfilename), dimensions);

    %%
    %{
    GROUND STATE SIMULATION TO GET STARTING FUNCTION PHI_0 FOR DYNAMIC SIMULATION
    %}

    %% Simulation methods

    Computation = 'Ground';
    Ncomponents = 1;
    Type = 'BESP';
    Deltat = 0.5;
    Stop_time = [];
    Stop_crit = {'MaxNorm', gammas(1)*1e-4};
    Max_iter = 1e3;

    Method = Method_Var3d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit, Max_iter);

    % temp save
    % Saving workspace with relevant data for fitting
    Method_ground = Method;
    save(info.get_workspace_path('fittingdata'), 'Method_ground');
    clear Method_ground;

    %% Geometry3D

% 	xmin = -2;
% 	xmax = 2;
% 	ymin= -1;
% 	ymax = 1;
% 	Nx = 2^8 + 1;
% 	Ny = 2^8 + 1;
    
    xmin = -xlimit;
    xmax = xlimit;
    ymin = -ylimit;
    ymax = ylimit;
    zmin = -zlimit;
    zmax = zlimit;
    Nx = xparticles;
    Ny = yparticles;
    Nz = zparticles;

    Geometry3D = Geometry3D_Var3d(xmin, xmax, ymin, ymax, zmin, zmax, Nx, Ny, Nz);

    %% Physics3D
    
    gx = gammas(1); gy = gammas(2); gz = gammas(3);

%     Delta = 0.5;
    Delta = delta;
    Beta = 4*pi*S;
    Omega = 0;
    Physics3D = Physics3D_Var3d(Method, Delta, Beta, Omega);
    %Physics3D = Potential_Var3d(Method, Physics3D); % std quadratic potential
        %Physics3D = Potential_Var3d(Method, Physics3D, @(X,Y,Z) quadratic_potential3d(1, 5, 1, X, Y, Z));

    %gx = 1; gy = 1; gz = 5;
    Physics3D = Potential_Var3d(Method, Physics3D, @(X,Y,Z) quadratic_potential3d(gx, gy, gz, X, Y, Z));
    Physics3D = Nonlinearity_Var3d(Method, Physics3D); % std cubic nonlinearity
    %Physics3D = Gradientx_Var3d(Method, Physics3D, @(x,y) -1i*Omega*y);
    %Physics3D = Gradienty_Var3d(Method, Physics3D, @(x,y) 1i*Omega*x);

    %% Defining a starting function Phi_0

    InitialData_choice = 1; % Gaussian initial data
    % w = (1 + Beta / (2*pi) )^(1/4); % interaction strength, w>1 for interactions; w=1 no interactions
    % s = 4; % size of the condensate for the trial function; s=1 is the expected result
    
    % find w3d
    exp = ExpPhis(chi, [gx,gy,gz], [0,0,0], Geometry3D);
    W = exp.getW();
    w = W.d3;
    clear exp W
    
    sigma_with_delta = struct();
    sigma_with_delta.x = w * (Delta/2)^(0.25) / sqrt(gx); sigma_with_delta.x = sigma_with_delta.x(sigma_with_delta.x>0);
    sigma_with_delta.y = w * (Delta/2)^(0.25) / sqrt(gy); sigma_with_delta.y = sigma_with_delta.y(sigma_with_delta.y>0);
    sigma_with_delta.z = w * (Delta/2)^(0.25) / sqrt(gz); sigma_with_delta.z = sigma_with_delta.z(sigma_with_delta.z>0);
    
    X0 = 0;
    Y0 = 0;
    Z0 = 0;
    gamma_x = 1;
    gamma_y = 1;
    gamma_z = 1;
    %gamma_x = 1 / (s*w^2);
    %gamma_y = 1 / (s*w^2);

    Phi_0 = InitialData_Var3d(Method, Geometry3D, Physics3D, InitialData_choice, X0, Y0, Z0, gamma_x, gamma_y, gamma_z);

    %% version mgmt
    curdir = strsplit(pwd, '/');
    curdir = strsplit(curdir{end}, '\');
    curdir = curdir{end};
    
    %% Printing interaction strength

    info.add_info_separator();
    info.add_custom_info('Folder: \t %s \n', curdir); % print current folder
    info.add_custom_info('Branch: \t %s \n', branch); % print current folder
    info.add_info_separator();
    info.add_custom_info('S \t=\t %f \n', S); % print interaction strength
    info.add_custom_info('Beta \t=\t %f \n', Beta); % print interaction parameter Beta
    info.add_custom_info('w \t=\t %f \n', w); % print Gaussian parameter w
    info.add_custom_info('sigmas \t=\t [%.4f,%.4f,%.4f] \n', ...
        sigma_with_delta.x,sigma_with_delta.y,sigma_with_delta.z); % print Gaussian parameter sigma_xyz with Delta
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
    Stop_crit = {'MaxNorm', 1e-4};

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

    %% Save PhiData structures for fitting w/o the whole workspace

    phi_dyn = PhiData(Phi, Geometry3D);
    phi_ground = PhiData(Phi_1, Geometry3D);
    phi_input = PhiData(Phi_0, Geometry3D);

    % Saving workspace with relevant data for fitting
    save(info.get_workspace_path('fittingdata'), 'gx', 'gy', 'gz', ... 
        'phi_dyn', 'phi_ground', 'phi_input', 'info', ... % necessary data
        'S', 'w', 'Beta', 'Delta', 'sigma_with_delta', 'Method', ... % additional data
        '-append'); % to not overwrite Method_ground

    simulation_finished = 'yes, but no plots yet!';
    save(info.get_workspace_path('fittingdata'),'simulation_finished', '-append')
    
    %% Some new output

    path = info.fulldir;
    wspath = info.get_workspace_path('fittingdata');
    datestr = info.creationTimeString;

    clearvars -except path wspath datestr dimensions

    load(wspath)

    PHIS = make_1d_cutouts_from_phi(phi_dyn.phisq{1});
    phix = PHIS{1}; phiy = PHIS{2}; phiz = PHIS{3};
    geom = phi_dyn.return_geometry();
    xx = geom.X(1,:,1);
    dx = geom.dx;
    %Redge = (15*S)^(1/5);
    area = dx*sum(phix);
    phixN = phix/area;

    save(wspath, 'PHIS', 'phix', 'phiy', 'phiz', 'geom', 'xx', 'dx', 'datestr', '-append')
    %save(wspath, 'Redge', 'area', 'phixN', '-append')
    save(wspath, 'area', 'phixN', '-append')

    expphis = ExpPhis(S, [gx,gy,gz], [0,0,0], geom);
    [g1d,wg1d] = expphis.Gaussian1d();
    [g2d,wg2d] = expphis.Gaussian2d();
    [g3d,wg3d] = expphis.Gaussian3d();
    [tf3d,redge3d] = expphis.TF3d();
    [tf2d,redge2d] = expphis.TF2d();
    [tf1d,redge1d] = expphis.TF1d();
    tfPHIS = make_1d_cutouts_from_phi(tf3d.phisq{1}); 
    tfphix = tfPHIS{1}; tfphiy = tfPHIS{2}; tfphiz = tfPHIS{3};
    tfphixN = tfphix / (dx*sum(tfphix));
    
    save(wspath,'tfphix','tfphixN','tfphiy', 'tfphiz', '-append')
    save(wspath,'expphis','g1d','g2d', 'g3d', 'tf1d','tf2d', 'tf3d', ...
        'redge1d', 'redge2d', 'redge3d', '-append')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yy = geom.Y(:,1,1)';
    zz = geom.Z(1,1,:); zz = zz(:);
    dy = geom.dy;
    dz = geom.dz;
    tfphiyN = tfphiy / (dy*sum(tfphiy));
    tfphizN = tfphiz / (dz*sum(tfphiz));
    gPHIS = make_1d_cutouts_from_phi(g3d.phisq{1}); gphix = gPHIS{1}; gphiy = gPHIS{2}; gphiz = gPHIS{3};
    gphixN = gphix / (dx*sum(gphix));
    gphiyN = gphiy / (dy*sum(gphiy));
    gphizN = gphiz / (dz*sum(gphiz));
    phiyN = phiy / (dy * sum(phiy));
    phizN = phiz / (dz * sum(phiz));

    save(wspath,'expphis','gphix','gphixN','gphiy','gphiyN', 'gphiz', 'gphizN', ...
        'tfphiyN', 'tfphizN', '-append')
    save(wspath,'yy','dy', 'zz', 'dz', 'phiyN', 'phizN', '-append')
    save(wspath, 'datestr', 'path', 'wspath', '-append')

    close all;
    %%%%%%%%%%%%%%%%%%%%% % set x boundaries
    if S >= 0.01
        nplots = 1;
        xlim1 = 3;
    elseif S < 0.01
        nplots = 3;
        xlim1 = 3;
        xlim2 = 1;
        xlim3 = 0.1;
        xlim4 = 0.2;
    end

    Rx = redge3d(1);
    Ry = redge3d(2);
    Rz = redge3d(3);
    
    %%%%%%%%%%%%%%%%%%%%%
    ylimit_x = max(max(max(gphixN),max(tfphixN)),max(phixN))*1.1;
    plot_1d_graph(xx,phixN,gphixN,tfphixN,xlim1,ylimit_x,'x',Rx,S,w,dx,datestr,Delta,dimensions)
    
    figname = ['compare-to-gaussTF-X' sprintf('_-%d;%d',xlim1,xlim1)];
    info.save_figure(1,'analyze',figname,path,'.fig')
    info.save_figure(1,'analyze',['compare-to-gaussTF-X' sprintf('_-%d;%d',xlim1,xlim1)],path,'.png')

    close all;
    %%%%%%%%%%%%%%%%%%%%%
    ylimit_y = max(max(max(gphiyN),max(tfphiyN)),max(phiyN))*1.1;
    plot_1d_graph(yy,phiyN,gphiyN,tfphiyN,xlim1,ylimit_y,'y',Ry,S,w,dy,datestr,Delta,dimensions)

    info.save_figure(1,'analyze',['compare-to-gaussTF-Y' sprintf('_-%d;%d',xlim1,xlim1)],path,'.fig')
    info.save_figure(1,'analyze',['compare-to-gaussTF-Y' sprintf('_-%d;%d',xlim1,xlim1)],path,'.png')
    
    close all;
    %%%%%%%%%%%%%%%%%%%%%
    ylimit_z = max(max(max(gphizN),max(tfphizN)),max(phizN))*1.1;
    plot_1d_graph(zz,phizN,gphizN,tfphizN,xlim1,ylimit_z,'z',Rz,S,w,dz,datestr,Delta,dimensions)

    info.save_figure(1,'analyze',['compare-to-gaussTF-Z' sprintf('_-%d;%d',xlim1,xlim1)],path,'.fig')
    info.save_figure(1,'analyze',['compare-to-gaussTF-Z' sprintf('_-%d;%d',xlim1,xlim1)],path,'.png')

    save(wspath, 'ylimit_x', 'ylimit_y', 'ylimit_z', '-append')
    
    %%%%%%%%%%%%%%%%%%%%%
    if nplots > 1
        close all;
        %%%%%%%%%%%%%%%%%%%%% % graph(-1,1)
        plot_1d_graph(xx,phixN,gphixN,tfphixN,xlim2,ylimit_x,'x',Rx,S,w,dx,datestr,Delta,dimensions)
        info.save_figure(1,'analyze',['compare-to-gaussTF-X' sprintf('_-%d;%d',xlim2,xlim2)],path,'.fig')
        info.save_figure(1,'analyze',['compare-to-gaussTF-X' sprintf('_-%d;%d',xlim2,xlim2)],path,'.png')

        close all;
        plot_1d_graph(yy,phiyN,gphiyN,tfphiyN,xlim2,ylimit_y,'y',Ry,S,w,dy,datestr,Delta,dimensions)
        info.save_figure(1,'analyze',['compare-to-gaussTF-Y' sprintf('_-%d;%d',xlim2,xlim2)],path,'.fig')
        info.save_figure(1,'analyze',['compare-to-gaussTF-Y' sprintf('_-%d;%d',xlim2,xlim2)],path,'.png')
        
        close all;
        plot_1d_graph(zz,phizN,gphizN,tfphizN,xlim2,ylimit_z,'z',Rz,S,w,dz,datestr,Delta,dimensions)
        info.save_figure(1,'analyze',['compare-to-gaussTF-Z' sprintf('_-%d;%d',xlim2,xlim2)],path,'.fig')
        info.save_figure(1,'analyze',['compare-to-gaussTF-Z' sprintf('_-%d;%d',xlim2,xlim2)],path,'.png')

        close all;
        %%%%%%%%%%%%%%%%%%%%% % graph(-0.1,0.1)
        plot_1d_graph(xx,phixN,gphixN,tfphixN,xlim3,ylimit_x,'x',Rx,S,w,dx,datestr,Delta,dimensions)
        info.save_figure(1,'analyze',['compare-to-gaussTF-X' sprintf('_-%d;%d',xlim3,xlim3)],path,'.fig')
        info.save_figure(1,'analyze',['compare-to-gaussTF-X' sprintf('_-%d;%d',xlim3,xlim3)],path,'.png')

        close all;
        plot_1d_graph(yy,phiyN,gphiyN,tfphiyN,xlim3,ylimit_y,'y',Ry,S,w,dy,datestr,Delta,dimensions)
        info.save_figure(1,'analyze',['compare-to-gaussTF-Y' sprintf('_-%d;%d',xlim3,xlim3)],path,'.fig')
        info.save_figure(1,'analyze',['compare-to-gaussTF-Y' sprintf('_-%d;%d',xlim3,xlim3)],path,'.png')
        
        close all;
        plot_1d_graph(zz,phizN,gphizN,tfphizN,xlim3,ylimit_z,'z',Rz,S,w,dz,datestr,Delta,dimensions)
        info.save_figure(1,'analyze',['compare-to-gaussTF-Z' sprintf('_-%d;%d',xlim3,xlim3)],path,'.fig')
        info.save_figure(1,'analyze',['compare-to-gaussTF-Z' sprintf('_-%d;%d',xlim3,xlim3)],path,'.png')
        
        close all;
        %%%%%%%%%%%%%%%%%%%%% % graph(-0.2,0.2)
        plot_1d_graph(xx,phixN,gphixN,tfphixN,xlim4,ylimit_x,'x',Rx,S,w,dx,datestr,Delta,dimensions)
        info.save_figure(1,'analyze',['compare-to-gaussTF-X' sprintf('_-%d;%d',xlim4,xlim4)],path,'.fig')
        info.save_figure(1,'analyze',['compare-to-gaussTF-X' sprintf('_-%d;%d',xlim4,xlim4)],path,'.png')

        close all;
        plot_1d_graph(yy,phiyN,gphiyN,tfphiyN,xlim4,ylimit_y,'y',Ry,S,w,dy,datestr,Delta,dimensions)
        info.save_figure(1,'analyze',['compare-to-gaussTF-Y' sprintf('_-%d;%d',xlim4,xlim4)],path,'.fig')
        info.save_figure(1,'analyze',['compare-to-gaussTF-Y' sprintf('_-%d;%d',xlim4,xlim4)],path,'.png')
        
        close all;
        plot_1d_graph(zz,phizN,gphizN,tfphizN,xlim4,ylimit_z,'z',Rz,S,w,dz,datestr,Delta,dimensions)
        
        info.save_figure(1,'analyze',['compare-to-gaussTF-Z' sprintf('_-%d;%d',xlim4,xlim4)],path,'.fig')
        info.save_figure(1,'analyze',['compare-to-gaussTF-Z' sprintf('_-%d;%d',xlim4,xlim4)],path,'.png')
    end

    simulation_finished = 'yes!';
    save(wspath,'simulation_finished', '-append')

    %% end
end
