function updatePlot(hObject, event, data, hplot)
%SETPLOTTIME Summary of this function goes here
%   Detailed explanation goes here
    
    if length(data) == 3
        Ldata = length(data1);
        flag_multiple = 1;
    else
        Ldata = length(data);
        flag_multiple = 0;
    end
    
    % Get the slider value and select the data
    sliderValue = round(get(hObject, 'Value'));
    if flag_multiple
        curData = cell(1, Ldata);
        for iplot = 1:Ldata
            curData{iplot} = data{iplot}{sliderValue};
        end
    else
        curData = data{sliderValue};
    end
    
    % Update the plot with the data
    if flag_multiple
        set(hplot, 'ydata', curData{1})
        hold on;
        set(hplot, 'ydata', curData{2})
        set(hplot, 'ydata', curData{3})
    else
        set(hplot, 'ydata', curData);   
    end
    
    % Determine limits
    upperlim = -Inf;
    lowerlim = Inf;
    for i = 1:Ldata
        % for long solution vectors we can disregard the first few in the
        % setting of the y-axis limits.
        if Ldata > 50 && i < 20
            continue;
        end
        upperlim = max(upperlim, max(data{i}));
        lowerlim = min(lowerlim, min(data{i}));
    end  

    % Picking the right limits (for the first 20 datasets)
    maxlim = max(upperlim, max(curData));
    minlim = min(lowerlim, min(curData));
    
    if abs(maxlim) < 10^(-10)
        maxlim = maxlim * 2/1.1;
    end
    if abs(minlim) < 10^(-10)
        minlim = minlim * 2/1.1;
    end
    
    % Scale axes: if one is dominant make sure to get a nice axis range
    if abs(maxlim/minlim) < 0.01
        maxlim = abs(minlim)/100;
    elseif abs(minlim/maxlim) < 0.01
        minlim = -abs(maxlim)/100;
    end
    
    % Scale axes so that there is space between data line & figure edges
    ylim([minlim*1.1 maxlim*1.1]);

    % Draw the plot!
    drawnow;
end

