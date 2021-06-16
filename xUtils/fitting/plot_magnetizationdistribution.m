% Plot population distributions of different mF in one figure on one axis

function [] = plot_magnetizationdistribution(geometry, Phi, info, direction)

    flag = []; % to specify whether direction was explicitly chosen
    if ~exist('direction','var') || isempty(direction)
        direction = 'x';
    else
        flag = 1; %direction explicitly given
    end
    
    close all;
    datestring = info.creationTimeString;
    
    % Get phi
    x = geometry.X;
    y = geometry.Y;
    z = geometry.Z;
    
    ly = floor(size(x,1)/2);
    lx = floor(size(y,2)/2);
    lz = floor(size(z,3)/2);
    
    M = abs(Phi{1}).^2 - abs(Phi{3}).^2;
    
    Mx = M(ly, :, lz);
    My = M(:, lx, lz);
    Mz = M(ly, lx, :); Mz = Mz(:)';
    
    % Find appropriate arrays
    if strcmp(direction, 'x')
        x = x(ly, :, lz);
        M = Mx;
        xax = 'x (centered in y,z) (a_{ho})';
    elseif strcmp(direction, 'y')
        y = y(:, lx, lz);
        M = My;
        x = y;
        xax = 'y (centered in x,z) (a_{ho})';
    elseif strcmp(direction, 'z')
        z = z(ly, lx, :); z = z(:)';
        M = Mz;
        x = z;
        xax = 'z (centered in x,y) (a_{ho})';
    else
        error('Something went wrong: direction is not recognized.')
    end
    
    limity = max(M);
    if limity < 10^(-10)
        limity = limity * 2/1.1;
    end
    
%     evomarker = floor(length(x)/30);
    evomarker = 1;
    % Creating figure
    plot(x, M, '--d', 'MarkerIndices', 1:evomarker:length(M), ...
        'LineWidth', 1.2, 'MarkerSize', 6)
    
    
%     plot(x, phi{1}, '--d', 'MarkerIndices', 1:evomarker:length(phi{1}), ...
%         'LineWidth', 1.2, 'MarkerSize', 6, 'Color', [0 0.6 0.1])
%     plot(x, phi{2}, '-.o', 'MarkerIndices', 1:evomarker:length(phi{2}), ...
%         'LineWidth', 1.2, 'MarkerSize', 6, 'Color', [0 0.2 0.2])
%     plot(x, phi{3}, '-.^', 'MarkerIndices', ceil(evomarker/2):evomarker:length(phi{3}), ...
%         'LineWidth', 1.2, 'MarkerSize', 6, 'Color', [0.8 0 0])
    
    ylim([-limity/100 limity*1.1]);
    
    % Add axes labels and figure text
    xlabel(xax); 
    ylabel('m = M/N');
    title('Magnetization distribution |\psi_+|^2-|\psi_-|^2');
    
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
    
    savename = 'Magnetization distribution';
    if flag
        savename = [savename '_' direction];
    end
    % Save figure
    info.save_figure(1, savename, '')
    hold off  
    
end