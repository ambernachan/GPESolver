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
    
    % Create the slider control
    theRange = length(solution) - 1;
    h = uicontrol('style','slider','units','pixel','position',[20 20 300 20],...
                  'min',1,'max',length(solution),'val',1,...
                  'sliderstep',[1/theRange 10/theRange]);
    
    % Create the plot
    hplot = plot(x, M{1}, '--d', 'MarkerIndices', 1:evomarker:L_yarray, ...
            'LineWidth', 1.2, 'MarkerSize', 6);
        
    % Register slider callback
    addlistener(h, 'ContinuousValueChange', @(hObject, event) updatePlot(hObject, event, M, hplot));

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

    % Save figure
    if flag % if the horizontal axis was explicitly chosen, put in savename
        savename = [savename '_' direction];
    end
    
    info.save_figure(1, savename, '')
    hold off
    
end