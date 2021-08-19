function [figHandles] = timeslider_slicer(geometry, solution, info, componentsPlotted)
    
    if ~exist('componentsPlotted', 'var')
        componentsPlotted = 1:length(solution{1});
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
    
    rho = cell(size(solution));
    for time = 1:length(solution)
        for component = 1:length(solution{1})
            rho{time}{component} = abs(solution{time}{component}).^2;
        end
    end
    
    for idx = 1:length(componentsPlotted)
        % create a time array of the component data
        datas = cell(1, length(rho));
        for time = 1:length(rho)
            datas{time} = rho{time}{idx};
        end
        
        % Create the slider and plot lines
        figNumber = idx;
        sliderHandle = createTimeSliderPlot(length(solution), figNumber);
        figHandle = gcf;
        figureName = ['abs(Phi_' str_mF{idx} ')^2'];
        figHandle.Name = figureName;
        lineHandle = slice(geometry.X, geometry.Y, geometry.Z, datas{1}, 0, 0, 0);
        
        % Register slider callback
        addlistener(sliderHandle, 'ContinuousValueChange', @(hObject, event) updateSlice(hObject, event, datas, lineHandle, geometry));
    %   
        colormap('jet')
        colorbar
        xlabel('x')
        ylabel('y')
        zlabel('z')
        str_title = sprintf('%s', str_mF{idx});
        title(['|\phi_{m_F=' str_title '}|^2'])

        % interpolate the colormap index value across the face of the plot,
        % colors become smooth
        shading(gca,'interp')
        
        % add simulation info and atom info
        datestring = info.creationTimeString;
        annotation('textbox', [0, 0.05, 0, 0], 'string', sprintf('%s', datestring))
        
        atom_str = getAtomStr(info.params.atom);

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