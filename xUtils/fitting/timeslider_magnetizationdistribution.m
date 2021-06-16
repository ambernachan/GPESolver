% Plot population distributions in one figure on one axis w/ time slider

function [] = timeslider_magnetizationdistribution(geometry, solution, info, direction)

    close all;
    datestring = info.creationTimeString;
    
    flag = []; % to specify whether direction was explicitly chosen
    if ~exist('direction','var') || isempty(direction)
        direction = 'x';
    else
        flag = 1; %direction explicitly given
    end

    % 3d spatial meshgrids
    x = geometry.X;
    y = geometry.Y;
    z = geometry.Z;
    
    % define midpoints for each of the x-arrays
    ly = floor(size(x,1)/2);
    lx = floor(size(y,2)/2);
    lz = floor(size(z,3)/2);
    
    % make the grid space axes 1d
    x = x(ly, :, lz);
    y = y(:, lx, lz);
    z = z(ly, lx, :); z = z(:)';
    
    % Find appropriate arrays: phi{time} 1d arrays
    for i = 1:length(solution)
        if strcmp(direction, 'x')
            Mx{i} = abs(solution{i}{1}(ly, :, lz)).^2 ...
                - abs(solution{i}{3}(ly, :, lz)).^2;
        elseif strcmp(direction, 'y')
            My{i} = abs(solution{i}{1}(:, lx, lz)).^2 ...
                - abs(solution{i}{3}(:, lx, lz)).^2;
        elseif strcmp(direction, 'z')
            Mz{i} = abs(solution{i}{1}(ly, lx, :)).^2 ...
                - abs(solution{i}{3}(ly, lx, :)).^2;
            Mz{i} = Mz{i}(:)';
        else
            error('Something went wrong: direction is not recognized.')
        end
    end
    
    if strcmp(direction, 'x')
        xax = 'x (centered in y,z) (a_{ho})';
        M = Mx;
    elseif strcmp(direction, 'y')
        xax = 'y (centered in x,z) (a_{ho})';
        x = y;
        M = My;
    elseif strcmp(direction, 'z')
        xax = 'z (centered in x,y) (a_{ho})';
        x = z;
        M = Mz;
    end
    
    limity = -Inf;
    lowerlim = Inf;
    for i = 1:length(solution)
        % for long solution vectors we can disregard the first few in the
        % setting of the y-axis limits.
        if length(solution) > 50 && i < 20
            continue;
        end
        limity = max(limity, max(M{i}));
        lowerlim = min(lowerlim, min(M{i}));
    end
    
    % Creating figure
    evomarker = floor(length(x)/20);
    L_yarray = length(M{1});
%     time = 1;
for time = 1:length(x)
    
    plot(x, M{time}, '--d', 'MarkerIndices', 1:evomarker:L_yarray, ...
        'LineWidth', 1.2, 'MarkerSize', 6)
    
    maxlim = max(limity, max(M{time}));
    minlim = min(lowerlim, min(M{time}));
    if maxlim < 10^(-10)
        maxlim = maxlim * 2/1.1;
    end
    if minlim < 0 && minlim > -10^(-10)
        minlim = minlim * 2/1.1;
        if abs(maxlim/minlim) < 0.01
            maxlim = abs(minlim)/100;
        end
    else
        minlim = -maxlim/100;
    end
    ylim([minlim maxlim*1.1]);
    
    % Add axes labels and figure text
    xlabel(xax); 
    ylabel('m = M/N');
    title('Magnetization distribution |\psi_+|^2-|\psi_-|^2');
    
%     lgd = legend('|\psi_+|', '|\psi_0|', '|\psi_-|');
    
    %add datestring to figure
    annotation('textbox', [0, 0.05, 0, 0], 'string', sprintf('%s', datestring))
    
    %add atom type text to figure
    wspath = info.get_workspace_path('groundstate');
    if ~exist(wspath, 'file')
        laptopstr = 'D:\surfdrive\UNI\MSc Project spinor condensates\';
        computerstr = 'C:\AMBER\surfdrive\UNI\MSc Project spinor condensates\';
        pathfromxOut = [];
        p = split(wspath, [string('/'), string('\')]);
        a = 0;
        for i=1:length(p)
            if strcmp(p(i), 'xOutputs')
                a = i;
            end
        end
        
        % In this case, the workspace is not defined on either PC/Laptop Amber,
        % you should save the figure manually.
        if a == 0
            sprintf('The atom type was not given because the workspace path was not recognized.')
            info.save_figure(1, 'Population distribution', '')
            hold off
            return;
        end
        
        p = p(a:end);
        for i=1:length(p)
            pathfromxOut = [pathfromxOut '\' char(p(i))];
        end
        
        if strcmp(wspath(1), 'C')
            wspath = [laptopstr pathfromxOut];
        elseif strcmp(wspath(1), 'D')
            wspath = [computerstr pathfromxOut];
        end
    end
    
    % Use the workspace path to load the atom mass, derive type of atom
    S = load(wspath, 'atom_mass'); atom_mass = S.atom_mass;
    atom_weight = atom_mass / getphysconst('amu');
    if atom_weight > 22 && atom_weight < 24
        atom_str = '^{23}Na';
    elseif atom_weight > 86 && atom_weight < 88
        atom_str = '^{87}Rb';
    else
        error('Please manually input atom type in plotting file.');
        atom_str = '';
    end
    
    % Add annotation about atom type to graph
    annotation('textbox', [0.15, 0.8, 0.1, 0.1], ...
        'string', sprintf('%s', atom_str), 'FitBoxToText', 'on')
    drawnow;
    hold off
end

    % Save figure
    savename = 'Magnetization distribution over time';
    if flag
        savename = [savename '_' direction];
    end
    
    info.save_figure(1, savename, '')
    hold off
    
end