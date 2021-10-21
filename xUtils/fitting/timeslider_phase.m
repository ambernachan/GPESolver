function [figHandles] = timeslider_phase(geometry, solution, info, componentsPlotted)
    
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
    
    angie = cell(size(solution));
    for time = 1:length(solution)
        for component = 1:length(solution{1})
            angie{time}{component} = angle(solution{time}{component});
        end
    end
    
    for idx = 1:length(componentsPlotted)
        % create a time array of the component data
        datas = cell(1, length(angie));
        for time = 1:length(angie)
            datas{time} = angie{time}{idx};
        end
        
        % Create the slider and plot lines
        figNumber = idx;
        sliderHandle = createTimeSliderPlot(length(solution), figNumber);
        figHandle = gcf;
        figureName = ['Phase of Phi{mF=' str_mF{idx} '}'];
        figHandle.Name = figureName;
        if isfield(geometry, 'Y')
            if isfield(geometry, 'Z')
                lineHandle = slice(geometry.X, geometry.Y, geometry.Z, datas{1}, 0, 0, 0);
            else
                error('2d plotting not yet implemented.')
            end
        else
            lineHandle = plot(geometry.X, datas{1});
        end
        
        % Register slider callback
        addlistener(sliderHandle, 'ContinuousValueChange', @(hObject, event) updateSlice(hObject, event, datas, lineHandle, geometry));
    %   
        colormap('jet')
        colorbar
        xlabel('x')
        ylabel('y')
        if isfield(geometry, 'Z')
            zlabel('z')
        end
        str_title = sprintf('%s', str_mF{idx});
        title(['Phase of |m_F= ' str_title '\rangle'])

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
        info.save_figure(idx, figureName, '')
        info.save_figure(idx, figureName, '', info.fulldir, '.png')
        hold off  

    end
    
    % get figure handles
    figHandles = flipud(findobj('Type', 'figure'));
    
end