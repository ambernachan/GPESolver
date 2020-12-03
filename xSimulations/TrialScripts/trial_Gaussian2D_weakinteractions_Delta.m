%{
...
%}

function [] = trial_Gaussian2D_weakinteractions_Delta(chi,xlimit,ylimit,xparticles,yparticles,delta)

%     clear;
%     chi = 0.01; xlimit = 3; xparticles = 2^5+1; ylimit = xlimit; yparticles = 2^5+1;
    close all;
    clearvars -except chi xlimit ylimit xparticles yparticles delta
    %clear;

    %% Determine interaction strength compared to kinetic energy

    S = chi;

    %% Saving info files

    dimensions = 2;
    info = Info(name_from_filename(mfilename), dimensions);

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
    Stop_crit = {'MaxNorm', 1e-4};
    Max_iter = 1e3;

    Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit, Max_iter);

    % temp save
    % Saving workspace with relevant data for fitting
    Method_ground = Method;
    save(info.get_workspace_path('fittingdata'), 'Method_ground');
    clear Method_ground;

    %% Geometry2D

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
    Nx = xparticles;
    Ny = yparticles;

    Geometry2D = Geometry2D_Var2d(xmin, xmax, ymin, ymax, Nx, Ny);

    %% Physics2D

%     Delta = 0.5;
    Delta = delta;
    Beta = 4*pi*S;
    Omega = 0;
    Physics2D = Physics2D_Var2d(Method, Delta, Beta, Omega);
    %Physics2D = Potential_Var2d(Method, Physics2D); % std quadratic potential
    %Physics2D = Potential_Var2d(Method, Physics2D, @(X,Y) quadratic_potential2d(1, 5, X, Y));
    Physics2D = Potential_Var2d(Method, Physics2D, @(X,Y) quadratic_potential2d(1, 1, X, Y));
    Physics2D = Nonlinearity_Var2d(Method, Physics2D); % std cubic nonlinearity
    %Physics2D = Gradientx_Var2d(Method, Physics2D, @(x,y) -1i*Omega*y);
    %Physics2D = Gradienty_Var2d(Method, Physics2D, @(x,y) 1i*Omega*x);

    %% Defining a starting function Phi_0

    InitialData_choice = 1; % Gaussian initial data
    w = (1 + Beta / (2*pi) )^(1/4); % interaction strength, w>1 for interactions; w=1 no interactions
    s = 4; % size of the condensate for the trial function; s=1 is the expected result
    X0 = 0;
    Y0 = 0;
    gamma_x = 1;
    gamma_y = 1;
    %gamma_x = 1 / (s*w^2);
    %gamma_y = 1 / (s*w^2);

    Phi_0 = InitialData_Var2d(Method, Geometry2D, Physics2D, InitialData_choice, X0, Y0, gamma_x, gamma_y);

    %% Printing interaction strength

    info.add_info_separator();
    info.add_custom_info('S \t=\t %f \n', S); % print interaction strength
    info.add_custom_info('Beta \t=\t %f \n', Beta); % print interaction parameter Beta
    info.add_custom_info('w \t=\t %f \n', w); % print Gaussian parameter w
    info.add_custom_info('Delta \t=\t %f \n', Delta); % print kinetic energy parameter Delta
    info.add_info_separator();

    %% Determining outputs

    Outputs = OutputsINI_Var2d(Method);

    %% Printing preliminary outputs

    Printing = 1;
    Evo = 15;
    Draw = 1;
    Print = Print_Var2d(Printing, Evo, Draw);

    %% RUN THE SIMULATION to find the ground state

    info.add_simulation_info(Geometry2D);
    [Phi_1, Outputs] = GPELab2d(Phi_0, Method, Geometry2D, Physics2D, Outputs, [], Print);

    %% Save the workspace & simulation info

    % save information about final iteration in info file
    info.add_result_info(Method, Outputs);
    % save workspace to workspace folder
    save(info.get_workspace_path('groundstate'));

    %% Draw & save solution

    close all;
    pause(2) % pauses the program for 2 seconds

    Draw_solution2d(Phi_0, Method, Geometry2D, Figure_Var2d());

    info.save_figure(1, 'initialdata', 'psi_sq',[],[]);
    info.save_figure(2, 'initialdata', 'angle',[],[]);

    Draw_solution2d(Phi_1, Method, Geometry2D, Figure_Var2d());

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
    Deltat = 1e-2;
    %Stop_time = 10;
    Max_iter = 200;
    Stop_crit = {'MaxNorm', 1e-4};

    %Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit, Max_iter);
    Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, [], Stop_crit, Max_iter);
    
    % Saving workspace with relevant data for fitting
    Method_dynamical = Method;
    save(info.get_workspace_path('fittingdata'), 'Method_dynamical', '-append');
    clear Method_dynamical;

    %% We keep Geometry2D as-is

    %% We keep Physics2D as-is

    %% Determining outputs

    Save_solution = 1;
    Outputs = OutputsINI_Var2d(Method, Save_solution);

    %% Printing preliminary outputs
    Printing = 1;
    Evo = 10;
    Draw = 1;
    Print = Print_Var2d(Printing, Evo, Draw);

    %% RUN THE SIMULATION to find the ground state

    [Phi, Outputs] = GPELab2d(Phi_1, Method, Geometry2D, Physics2D, Outputs, [], Print);

    %% Save the workspace & simulation info

    % save information about final simulation iteration in info file
    info.add_result_info(Method, Outputs);
    info.finish_info();
    save(info.get_workspace_path('dynamics'))

    %% Draw & save solution

    close all;
    pause(2) % pauses the program for 2 seconds

    Draw_solution2d(Phi, Method, Geometry2D, Figure_Var2d());

    info.save_figure(1, 'dynamics', 'psi_sq',[],[]);
    info.save_figure(2, 'dynamics', 'angle',[],[]);

    %% Save PhiData structures for fitting w/o the whole workspace

    phi_dyn = PhiData(Phi, Geometry2D);
    phi_ground = PhiData(Phi_1, Geometry2D);
    phi_input = PhiData(Phi_0, Geometry2D);

    % Saving workspace with relevant data for fitting
    save(info.get_workspace_path('fittingdata'), ... 
        'phi_dyn', 'phi_ground', 'phi_input', 'info', ... % necessary data
        'S', 'w', 'Beta', 'Delta', 'Method', ... % additional data
        '-append'); % to not overwrite Method_ground

    simulation_finished = 'yes, but no plots yet!';
    save(info.get_workspace_path('fittingdata'),'simulation_finished', '-append')
    
    %% Some new output

    path = info.fulldir;
    wspath = info.get_workspace_path('fittingdata');
    datestr = info.creationTimeString;

    clearvars -except path wspath datestr

    load(wspath)

    PHIS = make_1d_cutouts_from_phi(phi_dyn.phisq{1});
    phix = PHIS{1}; phiy = PHIS{2}';
    geom = phi_dyn.return_geometry();
    xx = geom.X(1,:);
    dx = geom.dx;
    Redge = (16*S)^(1/4);
    area = dx*sum(phix);
    phixN = phix/area;

    save(wspath, 'PHIS', 'phix', 'phiy', 'geom', 'xx', 'dx', 'datestr', '-append')
    save(wspath, 'Redge', 'area', 'phixN', '-append')

    expphis = ExpPhis(S, [1,1,1], [0,0,0], geom);
    [g1d,wg1d] = expphis.Gaussian1d();
    [g2d,wg2d] = expphis.Gaussian2d();
    [tf2d,redge2d] = expphis.TF2d();
    [tf1d,redge1d] = expphis.TF1d();
    tfPHIS = make_1d_cutouts_from_phi(tf2d.phisq{1}); tfphix = tfPHIS{1}; tfphiy = tfPHIS{2}';
    tfphixN = tfphix / (dx*sum(tfphix));
    
    save(wspath,'tfphix','tfphixN','tfphiy', '-append')
    save(wspath,'expphis','g1d','g2d','tf1d','tf2d', 'redge1d', 'redge2d','-append')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yy = geom.Y(:,1)';
    dy = geom.dy;
    tfphiyN = tfphiy / (dy*sum(tfphiy));
    gPHIS = make_1d_cutouts_from_phi(g2d.phisq{1}); gphix = gPHIS{1}; gphiy = gPHIS{2}';
    gphixN = gphix / (dx*sum(gphix));
    gphiyN = gphiy / (dy*sum(gphiy));
    phiyN = phiy / (dy * sum(phiy));

    save(wspath,'expphis','gphix','gphixN','gphiy','gphiyN','tfphiyN', '-append')
    save(wspath,'yy','dy', 'phiyN', '-append')
    save(wspath, 'datestr', 'path', '-append')

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

    %%%%%%%%%%%%%%%%%%%%%
    ylimit_x = max(max(max(gphixN),max(tfphixN)),max(phixN))*1.1;
    plot_1d_graph(xx,phixN,gphixN,tfphixN,xlim1,ylimit_x,'x',Redge,S,w,dx,datestr,Delta)
        
    figname = ['compare-to-gaussTF-X' sprintf('_-%d;%d',xlim1,xlim1)];
    info.save_figure(1,'analyze',figname,path,'.fig')
    info.save_figure(1,'analyze',['compare-to-gaussTF-X' sprintf('_-%d;%d',xlim1,xlim1)],path,'.png')

    close all;
    %%%%%%%%%%%%%%%%%%%%%
    ylimit_y = max(max(max(gphiyN),max(tfphiyN)),max(phiyN))*1.1;
    plot_1d_graph(yy,phiyN,gphiyN,tfphiyN,xlim1,ylimit_y,'y',Redge,S,w,dy,datestr,Delta)
    
    info.save_figure(1,'analyze',['compare-to-gaussTF-Y' sprintf('_-%d;%d',xlim1,xlim1)],path,'.fig')
    info.save_figure(1,'analyze',['compare-to-gaussTF-Y' sprintf('_-%d;%d',xlim1,xlim1)],path,'.png')

    save(wspath, 'ylimit_x', 'ylimit_y', '-append')
    
    %%%%%%%%%%%%%%%%%%%%%
    if nplots > 1
        close all;
        %%%%%%%%%%%%%%%%%%%%% % graph(-1,1)
        plot_1d_graph(xx,phixN,gphixN,tfphixN,xlim2,ylimit_x,'x',Redge,S,w,dx,datestr,Delta)
        info.save_figure(1,'analyze',['compare-to-gaussTF-X' sprintf('_-%d;%d',xlim2,xlim2)],path,'.fig')
        info.save_figure(1,'analyze',['compare-to-gaussTF-X' sprintf('_-%d;%d',xlim2,xlim2)],path,'.png')

        close all;
        plot_1d_graph(yy,phiyN,gphiyN,tfphiyN,xlim2,ylimit_y,'y',Redge,S,w,dy,datestr,Delta)
        info.save_figure(1,'analyze',['compare-to-gaussTF-Y' sprintf('_-%d;%d',xlim2,xlim2)],path,'.fig')
        info.save_figure(1,'analyze',['compare-to-gaussTF-Y' sprintf('_-%d;%d',xlim2,xlim2)],path,'.png')

        close all;
        %%%%%%%%%%%%%%%%%%%%% % graph(-0.1,0.1)
        plot_1d_graph(xx,phixN,gphixN,tfphixN,xlim3,ylimit_x,'x',Redge,S,w,dx,datestr,Delta)
        info.save_figure(1,'analyze',['compare-to-gaussTF-X' sprintf('_-%d;%d',xlim3,xlim3)],path,'.fig')
        info.save_figure(1,'analyze',['compare-to-gaussTF-X' sprintf('_-%d;%d',xlim3,xlim3)],path,'.png')

        close all;
        plot_1d_graph(yy,phiyN,gphiyN,tfphiyN,xlim3,ylimit_y,'y',Redge,S,w,dy,datestr,Delta)
        info.save_figure(1,'analyze',['compare-to-gaussTF-Y' sprintf('_-%d;%d',xlim3,xlim3)],path,'.fig')
        info.save_figure(1,'analyze',['compare-to-gaussTF-Y' sprintf('_-%d;%d',xlim3,xlim3)],path,'.png')
        
        close all;
        %%%%%%%%%%%%%%%%%%%%% % graph(-0.2,0.2)
        plot_1d_graph(xx,phixN,gphixN,tfphixN,xlim4,ylimit_x,'x',Redge,S,w,dx,datestr,Delta)
        info.save_figure(1,'analyze',['compare-to-gaussTF-X' sprintf('_-%d;%d',xlim4,xlim4)],path,'.fig')
        info.save_figure(1,'analyze',['compare-to-gaussTF-X' sprintf('_-%d;%d',xlim4,xlim4)],path,'.png')

        close all;
        plot_1d_graph(yy,phiyN,gphiyN,tfphiyN,xlim4,ylimit_y,'y',Redge,S,w,dy,datestr,Delta)
        info.save_figure(1,'analyze',['compare-to-gaussTF-Y' sprintf('_-%d;%d',xlim4,xlim4)],path,'.fig')
        info.save_figure(1,'analyze',['compare-to-gaussTF-Y' sprintf('_-%d;%d',xlim4,xlim4)],path,'.png')
    end

    simulation_finished = 'yes!';
    save(wspath,'simulation_finished', '-append')

    %% end
end
