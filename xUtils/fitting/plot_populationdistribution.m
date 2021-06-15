% Plot magnetism
function [] = plot_populationfractions(its, rho, info)

    close all;
    datestring = info.creationTimeString;
    x = 1:its;
    
    % Creating figure
    plot(x, rho{1}, '--')
    hold on
    plot(x, rho{2}, '-')
    plot(x, rho{3}, '-.')
    
    % Add axes labels and figure text
    xlabel('iterations'); 
    ylabel('|\phi_i|^2 / |\Psi|^2');
    title('Population fractions \rho_i');
    
    lgd = legend('\rho_+', '\rho_0', '\rho_-');
    
    %add datestring to figure
    annotation('textbox', [0, 0.05, 0, 0], 'string', sprintf('%s', datestring))
    
    %add atom type text to figure
    wspath = info.get_workspace_path('groundstate');
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
    
    annotation('textbox', [0.225, 0.2, 0.1, 0.1], ...
        'string', sprintf('%s', atom_str), 'FitBoxToText', 'on')
    
    % Save figure
    info.save_figure(1, 'Population fractions', '')
    hold off  
    
end