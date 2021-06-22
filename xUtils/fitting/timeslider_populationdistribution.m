% Plot population distributions in one figure on one axis w/ time slider

function [] = timeslider_populationdistribution(geometry, solution, info, direction)

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
            phiplus{i} = abs(solution{i}{1}(ly, :, lz));
            phizero{i} = abs(solution{i}{2}(ly, :, lz));
            phimin{i} = abs(solution{i}{3}(ly, :, lz));
        elseif strcmp(direction, 'y')
            phiplus{i} = abs(solution{i}{1}(:, lx, lz));
            phizero{i} = abs(solution{i}{2}(:, lx, lz));
            phimin{i} = abs(solution{i}{3}(:, lx, lz));
        elseif strcmp(direction, 'z')
            phiplus{i} = abs(solution{i}{1}(ly, lx, :)); phiplus{i} = phiplus{i}(:)';
            phizero{i} = abs(solution{i}{2}(ly, lx, :)); phizero{i} = phizero{i}(:)';
            phimin{i} = abs(solution{i}{3}(ly, lx, :)); phimin{i} = phimin{i}(:)';
        else
            error('Something went wrong: direction is not recognized.')
        end
    end
    
    if strcmp(direction, 'x')
        xax = 'x';
    elseif strcmp(direction, 'y')
        xax = 'y';
        x = y;
    elseif strcmp(direction, 'z')
        xax = 'z';
        x = z;
    end
    
    limity = -Inf;
    for i = 1:length(solution)
        % for long solution vectors we can disregard the first few in the
        % setting of the y-axis limits.
        if length(solution) > 50 && i < 20
            continue;
        end
        limity = max(limity, max(abs(phiplus{i})));
        limity = max(limity, max(abs(phizero{i})));
        limity = max(limity, max(abs(phimin{i})));
    end
    
    % Creating figure
    evomarker = floor(length(x)/20);
    L_yarray = length(phiplus{1});
%     time = 1;  

    % Create the slider and plot lines
    sliderHandle = createTimeSliderPlot(length(solution));
    hold on
    lh1 = plot(gca, x, phiplus{1}, '--d', 'MarkerIndices', 1:evomarker:L_yarray, ...
        'LineWidth', 1.2, 'MarkerSize', 6, 'Color', [0 0.6 0.1]);
    lh2 = plot(gca, x, phizero{1}, '-.o', 'MarkerIndices', 1:evomarker:L_yarray, ...
        'LineWidth', 1.2, 'MarkerSize', 6, 'Color', [0 0.2 0.2]);
    lh3 = plot(gca, x, phimin{1}, '-.^', 'MarkerIndices', ceil(evomarker/2):evomarker:L_yarray, ...
        'LineWidth', 1.2, 'MarkerSize', 6, 'Color', [0.8 0 0]);
    
    % Register slider callback
    addlistener(sliderHandle, 'ContinuousValueChange', @(hObject, event) updatePlot(hObject, event, {phiplus, phizero, phimin}, {lh1, lh2, lh3}));
%     t = timer('TimerFcn',@(x,y)disp('Hello World!'),'StartDelay',5);
    
    % Default axes: [left bottom width height] -> [0.1300 0.1100 0.7750 0.8150]
    set(gca, 'Position', [0.1300 0.1150 0.7750 0.8100]) 
    
    % Add axes labels and figure text
    xlabel(xax); 
    ylabel('|\phi_i|');
    title('Population distribution |\phi_i|');
    
    lgd = legend('|\psi_+|', '|\psi_0|', '|\psi_-|');
    
    savename = 'Population distribution over time';
    
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
    
    % Add datestring to figure
    annotation('textbox', [0, 0.05, 0, 0], 'string', sprintf('%s', datestring))

    % Add datestring to figure
    annotation('textbox', [0, 0.05, 0, 0], 'string', sprintf('%s', datestring))

    % if the horizontal axis was explicitly chosen, put in savename
    if flag 
        savename = [savename '_' direction];
    end
    
    % Save figure
    info.save_figure(1, savename, '')
    info.save_figure(1, savename, '', info.fulldir, '.png')
    hold off
    
end