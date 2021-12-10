function updatePlot(hObject, event, datas, lineHandles)
%SETPLOTTIME Summary of this function goes here
%   Detailed explanation goes here
    
    % Initialize upper and lower values of graph determination algorithm.
    maxval = -Inf;
    minval = Inf;
    
    % Get the slider value and select the data and line handle
    sliderValue = round(get(hObject, 'Value'));
    for i = 1:length(datas)
        curData = datas{i}{sliderValue};
        curLineHandle = lineHandles{i};
        % Update the line with the data
        set(curLineHandle, 'ydata', curData);
        
        % Determine graph limits by fancy algorithm
        for j = 1:length(datas{i})
            % for long solution vectors we can disregard the first few in the
            % setting of the y-axis limits.
            if length(datas{i}) > 50 && j < 20
                continue;
            end
            maxval = max(maxval, max(datas{i}{j}));
            minval = min(minval, min(datas{i}{j}));
        end  
    end
 
    % Picking the right limits (for the first 20 datasets)
    maxlim = max(maxval, max(curData));
    minlim = min(minval, min(curData));
    
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

