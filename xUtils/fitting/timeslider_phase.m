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
        lineHandle = slice(geometry.X, geometry.Y, geometry.Z, datas{1}, 0, 0, 0)
        
%         lh1 = plot(gca, x, phiplus{1}, '--d', 'MarkerIndices', 1:evomarker:L_yarray, ...
%             'LineWidth', 1.2, 'MarkerSize', 6, 'Color', [0 0.6 0.1]);
%         lh2 = plot(gca, x, phizero{1}, '-.o', 'MarkerIndices', 1:evomarker:L_yarray, ...
%             'LineWidth', 1.2, 'MarkerSize', 6, 'Color', [0 0.2 0.2]);
%         lh3 = plot(gca, x, phimin{1}, '-.^', 'MarkerIndices', ceil(evomarker/2):evomarker:L_yarray, ...
%             'LineWidth', 1.2, 'MarkerSize', 6, 'Color', [0.8 0 0]);

        % Register slider callback
        addlistener(sliderHandle, 'ContinuousValueChange', @(hObject, event) updateSlice(hObject, event, datas, lineHandle, geometry));
    %     t = timer('TimerFcn',@(x,y)disp('Hello World!'),'StartDelay',5);
       
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