% Plot population distributions of different mF in one figure on one axis

function [] = plot_populationdistribution(geometry, Phi, info, direction)

    if ~exist('direction','var') || isempty(direction)
        direction = 'x';
        dir = 1;
    else
        [direction, dir] = getDirection(direction);
    end
    xax = [direction ' (a_{ho})'];
    
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
    
    for n = 1:info.params.nComponents
        phi{n} = abs(Phi{n});
    end
    
    % Find appropriate arrays
    if dims == 3
        if strcmp(direction, 'x')
            x = X{1}(ly, :, lz);
            phix{1} = phi{1}(ly, :, lz);
            phix{2} = phi{2}(ly, :, lz);
            phix{3} = phi{3}(ly, :, lz);
            phi = phix;
        elseif strcmp(direction, 'y')
            y = X{2}(:, lx, lz);
            phiy{1} = phi{1}(:, lx, lz);
            phiy{2} = phi{2}(:, lx, lz);
            phiy{3} = phi{3}(:, lx, lz);
            phi = phiy;
            x = y;
        elseif strcmp(direction, 'z')
            z = X{3}(ly, lx, :); z = z(:)';
            phiz{1} = phi{1}(ly, lx, :); phiz{1} = phiz{1}(:)';
            phiz{2} = phi{2}(ly, lx, :); phiz{2} = phiz{2}(:)';
            phiz{3} = phi{3}(ly, lx, :); phiz{3} = phiz{3}(:)';
            phi = phiz;
            x = z;
        else
            error('Something went wrong: direction is not recognized.')
        end        
    elseif dims == 1
        x = X{1};
        xax = [direction ' (a_{ho})'];
        % phi doesn't need to change as it's already 1d
        [x, labl] = makexaxisinmeters(x, info);
        xax = [direction labl];
    else
        error('2-dimensional plotting not implemented yet.')
    end
       
    limity = max(max(max(phi{1}),max(phi{2})),max(phi{3}));
    
    % Creating figure
    evomarker = floor(length(x)/20);
    plot(x, phi{1}, '--d', 'MarkerIndices', 1:evomarker:length(phi{1}), ...
        'LineWidth', 1.2, 'MarkerSize', 6, 'Color', [0 0.6 0.1])
    hold on
    plot(x, phi{2}, '-.o', 'MarkerIndices', 1:evomarker:length(phi{2}), ...
        'LineWidth', 1.2, 'MarkerSize', 6, 'Color', [0 0.2 0.2])
    plot(x, phi{3}, '-.^', 'MarkerIndices', ceil(evomarker/2):evomarker:length(phi{3}), ...
        'LineWidth', 1.2, 'MarkerSize', 6, 'Color', [0.8 0 0])
    
    ylim([-limity/100 limity*1.1]);
    
    % Add axes labels and figure text
    xlabel(xax); 
    ylabel('|\phi_i|');
    title('Population distribution |\phi_i|');
    
    lgd = legend('|\psi_+|', '|\psi_0|', '|\psi_-|');
    
    savename = 'Population distribution';
    
    %add datestring to figure
    annotation('textbox', [0, 0.05, 0, 0], 'string', sprintf('%s', datestring))
    
    %add atom type text to figure
    atom_str = getAtomStr(info.params.atom);
    
    % Add annotation about atom type to graph
    annotation('textbox', [0.15, 0.8, 0.1, 0.1], ...
        'string', sprintf('%s', atom_str), 'FitBoxToText', 'on')
    
    % Save figure
    info.save_figure(1, savename, '')
    info.save_figure(1, savename, '', info.fulldir, '.png')
    hold off
    
end