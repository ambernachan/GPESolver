% Plot longitudinal magnetization
function [] = plot_longmagnetization(its, M, info, evo)

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
    
    % Creating figure
    evomarker = floor(length(x)/20);
    
    plot(x, M, '--d', 'MarkerIndices', 1:evomarker:length(M), ...
        'LineWidth', 1.2, 'MarkerSize', 6)
    
    % Add axes labels and figure text
    xlabel(labelx); 
    ylabel('m = M/N');
    title('Magnetization |\psi_+|^2-|\psi_-|^2');
    
%     lgd = legend('m');
    
    %add datestring to figure
    annotation('textbox', [0, 0.05, 0, 0], 'string', sprintf('%s', datestring))
    
    % gives wspath, which equals 0 when it hasn't been found.
    wspath = whichwspath(info); % gives wspath and parameter
    
    % Use the workspace path to load the atom mass, derive type of atom
    if wspath == 0 % meaning the wspath hasn't been found.
        info.save_figure(1, 'Magnetization', '')
        info.save_figure(1, 'Magnetization', '', info.fulldir, '.png')
        hold off
        return;
    end
    
    %add atom type text to figure
    atom_str = getAtomStr(info.params.atom);
    
    annotation('textbox', [0.725, 0.2, 0.1, 0.1], ...
        'string', sprintf('%s', atom_str), 'FitBoxToText', 'on')
    
    % Saving figure
    info.save_figure(1, 'Magnetization', '')
    info.save_figure(1, 'Magnetization', '', info.fulldir, '.png')
 
end