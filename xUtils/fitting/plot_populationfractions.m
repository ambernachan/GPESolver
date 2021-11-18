% Plot magnetism
function [] = plot_populationfractions(its, rho, info, evo, method)

    close all;
    datestring = info.creationTimeString;
    if ~exist('evo','var') || isempty(evo)
        evo = 1;
    end
    
    % Determine whether this was a ground state or dynamic simulation
    % if 'ground' exists in the path name; OR neither 'ground' nor 'dynamic' does
    if ~isempty(strfind(info.name, 'ground')) || ...
            (isempty(strfind(info.name, 'ground')) && isempty(strfind(info.name, 'dynamic')))
        labelx = 'iterations';
        x = evo:evo:its*evo;
    end
    % if 'ground' doesn't exist in the path name but 'dynamic' does
    if isempty(strfind(info.name, 'ground'))
        if ~isempty(strfind(info.name, 'dynamic'))
            labelx = 'time (\omega_{osc}^{-1})';
            if exist('method', 'var')
                dt = method.Deltat;
            else
                dt = 1;
            end
            x = evo*dt:evo*dt:its*evo*dt; % time in 1/w_ho
            [x, labelx] = makedynamictimeline(x, info);
        end
    end
    
%     x = evo:evo:its*evo;
    
    limity = max(max(max(rho{1}),max(rho{2})),max(rho{3}));
    lowery = min(min(min(rho{1}),min(rho{2})),min(rho{3}));
    
    % Creating figure
    evomarker = floor(length(x)/10);
    plot(x, rho{1}, '--d', 'MarkerIndices', 1:evomarker:length(rho{1}), ...
        'LineWidth', 1.2, 'MarkerSize', 6, 'Color', [0 0.6 0.1])
    hold on
    plot(x, rho{2}, '-.o', 'MarkerIndices', 1:evomarker:length(rho{2}), ...
        'LineWidth', 1.2, 'MarkerSize', 6, 'Color', [0 0.2 0.2])
    plot(x, rho{3}, '-.^', 'MarkerIndices', ceil(evomarker/2):evomarker:length(rho{3}), ...
        'LineWidth', 1.2, 'MarkerSize', 6, 'Color', [0.8 0 0])
    
    ylim([lowery-limity/100 limity*1.1]);
    
    % Add axes labels and figure text
    xlabel(labelx); 
    ylabel('|\phi_i|^2 / |\Psi|^2');
    title('Population fractions \rho_i');
    
    lgd = legend('\rho_+', '\rho_0', '\rho_-');
    
    %add datestring to figure
    annotation('textbox', [0, 0.05, 0, 0], 'string', sprintf('%s', datestring))
    
    savename = 'Population fractions';
    
    %add atom type text to figure
    atom_str = getAtomStr(info.params.atom);
    
%     annotation('textbox', [0.225, 0.2, 0.1, 0.1], ...
%         'string', sprintf('%s', atom_str), 'FitBoxToText', 'on')
    annotation('textbox', [0.4, 0.2, 0.1, 0.1], ...
        'string', sprintf('%s', atom_str), 'FitBoxToText', 'on', ...
        'BackgroundColor', [1, 1, 1], 'FaceAlpha', 0.8)
    
    % Save figure
    info.save_figure(1, savename, '')
    info.save_figure(1, savename, '', info.fulldir, '.png')
    hold off  
    
end