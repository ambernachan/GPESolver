function [figHandles] = plot_phase(geometry, Phi, info, componentsPlotted)
    
    if ~exist('componentsPlotted', 'var')
        componentsPlotted = 1:length(Phi);
    end
    
    for i = 1:length(componentsPlotted)
        if componentsPlotted(i) == 1
            str_mF{i} = '1';
        elseif componentsPlotted(i) == 2
            str_mF{i} = '0';
        elseif componentsPlotted(i) == 3
            str_mF{i} = '-1';
        end
    end
    
    for idx = 1:length(componentsPlotted)
        figHandle = figure(idx);
        figureName = ['Phase of Phi{mF=' str_mF{idx} '}'];
        figHandle.Name = figureName;
        slice(geometry.X, geometry.Y, geometry.Z, angle(Phi{idx}), 0, 0, 0)
        colormap('jet')
        colorbar
        xlabel('x')
        ylabel('y')
        zlabel('z')
        str_title = sprintf('%s', str_mF{idx});
%         title(['Phase of \phi_{mF= ' str_title '}'])
        title(['Phase of |m_F= ' str_title '\rangle'])

        % interpolate the colormap index value across the face of the plot, colors
        % become smooth
        shading(gca,'interp')
        
        % add simulation info and atom info
        datestring = info.creationTimeString;
        annotation('textbox', [0, 0.05, 0, 0], 'string', sprintf('%s', datestring))
        
        if ~isfield(info.params, 'atom')
            info.save_figure(1, figureName, '')
            info.save_figure(1, figureName, '', info.fulldir, '.png')
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

        annotation('textbox', [0.1375, 0.825, 0.1, 0.1], ...
            'string', sprintf('%s', atom_str), 'FitBoxToText', 'on')

        % Save figure
        info.save_figure(1, figureName, '')
        info.save_figure(1, figureName, '', info.fulldir, '.png')
        hold off  

    end
    
    % get figure handles
    figHandles = flipud(findobj('Type', 'figure'));
    
end