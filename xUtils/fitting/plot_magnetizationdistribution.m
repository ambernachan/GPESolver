% Plot population distributions of different mF in one figure on one axis

function [] = plot_magnetizationdistribution(geometry, Phi, info, direction)

    flag = []; % to specify whether direction was explicitly chosen
    if ~exist('direction','var') || isempty(direction)
        direction = 'x';
        dir = 1;
    else
        flag = 1; %direction explicitly given
        [direction, dir] = getDirection(direction);
    end
    xax = direction;
    
    dims = getDimensionality(geometry);
    directions = [{'X'}, {'Y'}, {'Z'}];
    
    close all;
    datestring = info.creationTimeString;
    
    % Get phi
    for d=1:dims
        X{d} = geometry.(directions{d});
    end
    
    if dims == 1
        lx = floor(numel(X{1})/2);
    else
        ly = floor(size(X{1},1)/2);
        lx = floor(size(X{2},2)/2);
        if dims == 3
            lz = floor(size(X{3},3)/2);
        end
    end
    
    M = abs(Phi{1}).^2 - abs(Phi{3}).^2;
    
    if dims > 1
        Mx = M(ly, :, lz);
        My = M(:, lx, lz);
        if dims > 2
            Mz = M(ly, lx, :); Mz = Mz(:)';
        end
    elseif dims == 1
        Mx = M;
    end
    
    % Sanity check for matching directions of plotting to dimensionality
    if dir == 3
        if dims < 3
            error('Dimensions required do not match direction specified.')
        end
    elseif dir == 2
        if dims < 2
            error('Dimensions required do not match direction specified.')
        end
    end
    if dims == 1
        if dir>1
            error('Dimensions required do not match direction specified.')
        end
    elseif dims == 2
        if dir>2
            error('Dimensions required do not match direction specified.')
        end
    end
    
    if dir == 1 % x
        x = X{1};
        % M = Mx already holds
    elseif dir == 2 % y
        if dims == 2
            error('2d plotting not yet supported for magnetizationdistribution')
        elseif dims == 3
            y = y(:, lx, lz);
            M = My;
            x = y;
            xax = 'y (centered in x,z) (a_{ho})';
        end
    elseif dir == 3 % z
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
%         ylim([-limity/100 limity*1.1]);
        ylim ([-1.1 1.1]);
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
    
    % Use the workspace path to load the atom mass, derive type of atom
    if wspath == 0 % meaning the wspath hasn't been found.
        info.save_figure(1, savename, '')
        info.save_figure(1, savename, '', info.fulldir, '.png')
        hold off
        return;
    end
    
    %add atom type text to figure
    atom_str = getAtomStr(info.params.atom);
    
    % Add annotation about atom type to graph
%     annotation('textbox', [0.15, 0.8, 0.1, 0.1], ...
%         'string', sprintf('%s', atom_str), 'FitBoxToText', 'on')
    annotation('textbox', [0.4, 0.2, 0.1, 0.1], ...
        'string', sprintf('%s', atom_str), 'FitBoxToText', 'on', ...
        'BackgroundColor', [1, 1, 1], 'FaceAlpha', 0.8)
    
    if flag
        savename = [savename '_' direction];
    end
    % Save figure
    info.save_figure(1, savename, '')
    info.save_figure(1, savename, '', info.fulldir, '.png')
    hold off  
    
end