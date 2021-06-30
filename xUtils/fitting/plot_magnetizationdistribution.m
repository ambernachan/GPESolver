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
    
    if limity < 0
        lowerlim = min(M);
        if abs(lowerlim) > abs(limity)
            ylim([lowerlim*1.1 2*limity*1.1]);
        else
            ylim([lowerlim*1.1 limity*1.1]);
        end
    else
        ylim([-limity/100 limity*1.1]);
    end
    
    % Add axes labels and figure text
    xlabel(xax); 
    ylabel('m = M/N');
    title('Magnetization distribution |\psi_+|^2-|\psi_-|^2');
    
    %add datestring to figure
    annotation('textbox', [0, 0.05, 0, 0], 'string', sprintf('%s', datestring))
    
    % gives wspath, which equals 0 when it hasn't been found.
    wspath = whichwspath(info); % gives wspath and parameter
    
    savename = 'Magnetization distribution';
    
    % add atom type text to figure
    % Use the workspace path to load the atom mass, derive type of atom
    if wspath == 0 % meaning the wspath hasn't been found.
        info.save_figure(1, savename, '')
        info.save_figure(1, savename, '', info.fulldir, '.png')
        hold off
        return;
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
    
    if flag
        savename = [savename '_' direction];
    end
    % Save figure
    info.save_figure(1, savename, '')
    info.save_figure(1, savename, '', info.fulldir, '.png')
    hold off  
    
end