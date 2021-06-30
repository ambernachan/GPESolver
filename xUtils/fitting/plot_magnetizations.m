% Plot magnetization, M = [Mx; My; Mz]

function [] = plot_magnetizations(its, M, info, evo)

    close all;
    datestring = info.creationTimeString;
    if ~exist('evo','var') || isempty(evo)
        evo = 1;
    end
    
    x = evo:evo:its*evo;
    
    % marker setting step size
    evomarker = floor(length(x)/20);
    % Creating figure
    plot(x, M(1, :), '--x', 'MarkerIndices', 1:evomarker:length(M(1)), ...
        'LineWidth', 1.2, 'MarkerSize', 6, 'Color', [0 0.6 0.1])
    hold on
    plot(x, M(2, :), '-.+', 'MarkerIndices', 1:evomarker:length(M(2)), ...
        'LineWidth', 1.2, 'MarkerSize', 6, 'Color', [0 0.2 0.2])
    plot(x, M(3, :), '-.o', 'MarkerIndices', ceil(evomarker/2):evomarker:length(M(3)), ...
        'LineWidth', 1.2, 'MarkerSize', 6, 'Color', [0.8 0 0])
    
    % Add axes labels and figure text
    xlabel('iterations'); 
    ylabel('m = M/N');
    title('Magnetization <Fx>, <Fy>, <Fz>');
    
    lgd = legend('<Fx>', '<Fy>', '<Fz>');
    
    %add datestring to figure
    annotation('textbox', [0, 0.05, 0, 0], 'string', sprintf('%s', datestring))
    
    savename = 'Magnetization in x,y,z';
    
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
    
%     annotation('textbox', [0.725, 0.2, 0.1, 0.1], ...
%         'string', sprintf('%s', atom_str), 'FitBoxToText', 'on')
    annotation('textbox', [0.225, 0.2, 0.1, 0.1], ...
        'string', sprintf('%s', atom_str), 'FitBoxToText', 'on')
    
    % Save figure
    info.save_figure(1, savename, '')
    info.save_figure(1, savename, '', info.fulldir, '.png')
    hold off
 
end