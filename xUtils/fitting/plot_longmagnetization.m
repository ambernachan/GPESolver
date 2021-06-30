% Plot longitudinal magnetization
function [] = plot_longmagnetization(its, M, info, evo)

    close all;
    datestring = info.creationTimeString;
    if ~exist('evo','var') || isempty(evo)
        evo = 1;
    end
    
    x = evo:evo:its*evo;
    
    % Creating figure
    evomarker = floor(length(x)/20);
    
    plot(x, M, '--d', 'MarkerIndices', 1:evomarker:length(M), ...
        'LineWidth', 1.2, 'MarkerSize', 6)
    
    % Add axes labels and figure text
    xlabel('iterations'); 
    ylabel('m = M/N');
    title('Magnetization |\psi_+|^2-|\psi_-|^2');
    
%     lgd = legend('m');
    
    %add datestring to figure
    annotation('textbox', [0, 0.05, 0, 0], 'string', sprintf('%s', datestring))
    
    % gives wspath, which equals 0 when it hasn't been found.
    wspath = whichwspath(info); % gives wspath and parameter
    
    %add atom type text to figure
    % Use the workspace path to load the atom mass, derive type of atom
    if wspath == 0 % meaning the wspath hasn't been found.
        info.save_figure(1, 'Magnetization', '')
        info.save_figure(1, 'Magnetization', '', info.fulldir, '.png')
        hold off
        return;
    end
    
    S = load(wspath, 'atom_mass'); atom_mass = S.atom_mass;
    atom_weight = atom_mass / getphysconst('amu');
    if atom_weight > 22 && atom_weight < 24
        atom_str = '^{23}Na';
    elseif atom_weight > 86 && atom_weight < 88
        atom_str = '^{87}Rb';
    else
        error('Please manually input atom type in plotting file.');
        atom_str = '';
    end
    
    annotation('textbox', [0.725, 0.2, 0.1, 0.1], ...
        'string', sprintf('%s', atom_str), 'FitBoxToText', 'on')
    
    % Saving figure
    info.save_figure(1, 'Magnetization', '')
    info.save_figure(1, 'Magnetization', '', info.fulldir, '.png')
 
end