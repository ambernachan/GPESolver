% Plot population distributions in one figure on one axis w/ time slider

function [] = timeslider_dirmagnetizationdistribution(geometry, solution, info, direction)

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
    
    Ncomponents = 3; method.Ncomponents = Ncomponents;
    % Find appropriate arrays: phi{time} 1d arrays
    for iter = 1:length(solution)
        Fdistr = findExpFdistr(method, geometry, solution{iter});
        for n = 1:Ncomponents
            xF{iter}{n} = sum(sum(Fdistr{n}, 1), 3);
            yF{iter}{n} = sum(sum(Fdistr{n}, 2), 3);
            zF{iter}{n} = sum(sum(Fdistr{n}, 1), 2);
        end
        
        for n = 1:Ncomponents
            phisq_x{iter}{n} = 0;
            phisq_y{iter}{n} = 0;
            phisq_z{iter}{n} = 0;
            
            for iy = 1:size(x,1)
                for iz = 1:size(z,3)
                    phisq_x{iter}{n} = phisq_x{iter}{n} + abs(solution{iter}{n}(iy, :, iz)).^2;
                end
            end
            for ix = 1:size(y,2)
                for iz = 1:size(z,3)
                    phisq_y{iter}{n} = phisq_y{iter}{n} + abs(solution{iter}{n}(:, ix, iz)).^2;
                end
            end
            for ix = 1:size(y,2)
                for iy = 1:size(x,1)
                    phisq_z{iter}{n} = phisq_z{iter}{n} + abs(solution{iter}{n}(iy, ix, :)).^2;
                end
            end
            phisq_z{iter}{n} = phisq_z{iter}{n}(:)';
        end

        if strcmp(direction, 'x')
            Fz{iter} = phisq_x{iter}{1} - phisq_x{iter}{3};
        elseif strcmp(direction, 'y')
            Fz{iter} = phisq_y{iter}{1} - phisq_y{iter}{3};
        elseif strcmp(direction, 'z')
            Fz{iter} = phisq_z{iter}{1} - phisq_z{iter}{3};
            Fz{iter} = Fz{iter}(:)';
        else
            error('Something went wrong: direction is not recognized.')
        end
    end
    
    % Create 1d plot data: ydata M set to appropriate direction & xaxis
    % title selected
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
     
    % Creating figure
    evomarker = floor(length(x)/20);
    L_yarray = length(M{1});
    
    % Create the slider and plot lines
    sliderHandle = createTimeSliderPlot(length(solution));
    lh = plot(gca, x, M{1}, '--d', 'MarkerIndices', 1:evomarker:L_yarray, ...
            'LineWidth', 1.2, 'MarkerSize', 6);
        
    % Register slider callback
    addlistener(sliderHandle, 'ContinuousValueChange', @(hObject, event) updatePlot(hObject, event, {M}, {lh}));

    % Default axes: [left bottom width height] -> [0.1300 0.1100 0.7750 0.8150]
    set(gca, 'Position', [0.1300 0.1150 0.7750 0.8100]) 
    
    % Add axes labels and figure text
    xlabel(xax); 
    ylabel('m = M/N');
    title('Magnetization distribution |\psi_+|^2-|\psi_-|^2');

    savename = 'Magnetization distribution over time';

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

    % if the horizontal axis was explicitly chosen, put in savename
    if flag 
        savename = [savename '_' direction];
    end
    
    % Save figure
    info.save_figure(1, savename, '')
    info.save_figure(1, savename, '', info.fulldir, '.png')
    hold off
    
end