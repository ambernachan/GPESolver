% Plot population distributions in one figure on one axis w/ time slider

function [] = timeslider_populationdistribution(geometry, solution, info, direction)

    close all;
    datestring = info.creationTimeString;
    
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

    % 3d spatial meshgrids
    for d=1:dims
        X{d} = geometry.(directions{d});
    end
    
    % define midpoints for each of the x-arrays
    if dims == 1
        lx = floor(numel(X{1})/2);
    else
        ly = floor(size(X{1},1)/2);
        lx = floor(size(X{2},2)/2);
        if dims == 3
            lz = floor(size(X{3},3)/2);
        end
    end

    % make the grid space axes 1d
    if dims == 3
        if strcmp(direction, 'x')
            x = X{1}(ly, :, lz);
        elseif strcmp(direction, 'y')
            y = X{2}(:, lx, lz);
            x = y;
        elseif strcmp(direction, 'z')
            z = X{3}(ly, lx, :); z = z(:)';
            x = z;
        else
            error('Something went wrong: direction is not recognized.')
        end        
    elseif dims == 1
        x = X{1};
    else
        error('2-dimensional plotting not implemented yet.')
    end
    
    % Find appropriate arrays: phi{time} 1d arrays
    if dims == 3
        for i = 1:length(solution)
            if strcmp(direction, 'x')
                phiplus{i} = abs(solution{i}{1}(ly, :, lz));
                phizero{i} = abs(solution{i}{2}(ly, :, lz));
                phimin{i} = abs(solution{i}{3}(ly, :, lz));
            elseif strcmp(direction, 'y')
                phiplus{i} = abs(solution{i}{1}(:, lx, lz));
                phizero{i} = abs(solution{i}{2}(:, lx, lz));
                phimin{i} = abs(solution{i}{3}(:, lx, lz));
                x = y;
            elseif strcmp(direction, 'z')
                phiplus{i} = abs(solution{i}{1}(ly, lx, :)); phiplus{i} = phiplus{i}(:)';
                phizero{i} = abs(solution{i}{2}(ly, lx, :)); phizero{i} = phizero{i}(:)';
                phimin{i} = abs(solution{i}{3}(ly, lx, :)); phimin{i} = phimin{i}(:)';
                x = z;
            else
                error('Something went wrong: direction is not recognized.')
            end
        end
    elseif dims == 1
        for i = 1:length(solution)
            phiplus{i} = abs(solution{i}{1});
            phizero{i} = abs(solution{i}{2});
            phimin{i} = abs(solution{i}{3});
        end
    elseif dims == 2
        error('2-dimensional plotting not implemented yet.')
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
    
    %add atom type text to figure
    atom_str = getAtomStr(info.params.atom);
    
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