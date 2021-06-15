% Plot magnetism
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
    elseif strcmp(direction, 'y')
        y = y(:, lx, lz);
        phiy{1} = phi{1}(:, lx, lz);
        phiy{2} = phi{2}(:, lx, lz);
        phiy{3} = phi{3}(:, lx, lz);
        phi = phiy;
        x = y;
    elseif strcmp(direction, 'z')
        z = z(ly, lx, :); z = z(:)';
        phiz{1} = phi{1}(ly, lx, :); phiz{1} = phiz{1}(:)';
        phiz{2} = phi{2}(ly, lx, :); phiz{2} = phiz{2}(:)';
        phiz{3} = phi{3}(ly, lx, :); phiz{3} = phiz{3}(:)';
        phi = phiz;
        x = z;
    else
        error('Something went wrong: direction is not recognized.')
    end
    
    % Creating figure
    plot(x, phi{1}, '--')
    hold on
    plot(x, phi{2}, '-')
    plot(x, phi{3}, '-.')
    
    % Add axes labels and figure text
    xlabel('x'); 
    ylabel('|\phi_i|');
    title('Population distribution |\phi_i|');
    
    lgd = legend('\psi_+', '\psi_0', '\psi_-');
    
    %add datestring to figure
    annotation('textbox', [0, 0.05, 0, 0], 'string', sprintf('%s', datestring))
    
    %add atom type text to figure
    wspath = info.get_workspace_path('groundstate');
    if ~exist(wspath, 'file')
        laptopstr = 'D:\surfdrive\UNI\MSc Project spinor condensates\';
        computerstr = 'C:\AMBER\surfdrive\UNI\MSc Project spinor condensates\';
        pathfromxOut = [];
        p = split(wspath, [string('/'), string('\')]);
        a = 0;
        for i=1:length(p)
            if strcmp(p(i), 'xOutputs')
                a = i;
            end
        end
        
        % In this case, the workspace is not defined on either PC/Laptop Amber,
        % you should save the figure manually.
        if a == 0
            sprintf('The atom type was not given because the workspace path was not recognized.')
            info.save_figure(1, 'Population distribution', '')
            hold off
            return;
        end
        
        p = p(a:end);
        for i=1:length(p)
            pathfromxOut = [pathfromxOut '\' char(p(i))];
        end
        
        if strcmp(wspath(1), 'C')
            wspath = [laptopstr pathfromxOut];
        elseif strcmp(wspath(1), 'D')
            wspath = [computerstr pathfromxOut];
        end
    end
    
    % Use the workspace path to load the atom mass, derive type of atom
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
    
    % Add annotation about atom type to graph
    annotation('textbox', [0.225, 0.75, 0.1, 0.1], ...
        'string', sprintf('%s', atom_str), 'FitBoxToText', 'on')
    
    % Save figure
    info.save_figure(1, 'Population distribution', '')
    hold off  
    
end