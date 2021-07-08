% Plot population distributions of different mF in one figure on one axis

function [] = plot_populationdistribution(geometry, Phi, info, direction)


    if ~exist('direction','var') || isempty(direction)
        direction = 'x';
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
    
    phi{1} = abs(Phi{1});
    phi{2} = abs(Phi{2});
    phi{3} = abs(Phi{3});
    
    % Find appropriate arrays
    if strcmp(direction, 'x')
        x = x(ly, :, lz);
        phix{1} = phi{1}(ly, :, lz);
        phix{2} = phi{2}(ly, :, lz);
        phix{3} = phi{3}(ly, :, lz);
        phi = phix;
        xax = 'x';
    elseif strcmp(direction, 'y')
        y = y(:, lx, lz);
        phiy{1} = phi{1}(:, lx, lz);
        phiy{2} = phi{2}(:, lx, lz);
        phiy{3} = phi{3}(:, lx, lz);
        phi = phiy;
        x = y;
        xax = 'y';
    elseif strcmp(direction, 'z')
        z = z(ly, lx, :); z = z(:)';
        phiz{1} = phi{1}(ly, lx, :); phiz{1} = phiz{1}(:)';
        phiz{2} = phi{2}(ly, lx, :); phiz{2} = phiz{2}(:)';
        phiz{3} = phi{3}(ly, lx, :); phiz{3} = phiz{3}(:)';
        phi = phiz;
        x = z;
        xax = 'z';
    else
        error('Something went wrong: direction is not recognized.')
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
    
    %add datestring to figure
    annotation('textbox', [0, 0.05, 0, 0], 'string', sprintf('%s', datestring))
    
    if ~isfield(info.params, 'atom')
        info.save_figure(1, savename, '')
        info.save_figure(1, savename, '', info.fulldir, '.png')
        sprintf('Warning: atom type was not specified')
        hold off
        return;
    end

    if strcmp(info.params.atom, 'Na')
        atom_str = '^{23}Na';
    elseif strcmp(info.params.atom, 'Rb')
        atom_str = '^{87}Rb';
    else
        atom_str = info.params.atom;
    end
    
    % Add annotation about atom type to graph
    annotation('textbox', [0.15, 0.8, 0.1, 0.1], ...
        'string', sprintf('%s', atom_str), 'FitBoxToText', 'on')
    
    % Save figure
    info.save_figure(1, savename, '')
    info.save_figure(1, savename, '', info.fulldir, '.png')
    hold off
    
end