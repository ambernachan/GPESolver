function updateSlice(hObject, event, datas, lineHandle, geometry)
%SETPLOTTIME Summary of this function goes here
%   Detailed explanation goes here
    
    % Get the slider value and select the data and line handle
    sliderValue = round(get(hObject, 'Value'));
    curData = datas{sliderValue};
    
    % create new (temporary) axes
    num = 99;
    newfig = figure(num);
    newax = axes(newfig);
    % plot the sliced angle into the temporary axes
    if isfield(geometry, 'Y')
        if isfield(geometry, 'Z')
            curSlicer = slice(newax, geometry.X,geometry.Y,geometry.Z, curData, 0,0,0);
        end
    else
        curSlicer = plot(newax, geometry.X, curData);
    end
    
    curLineHandle = lineHandle;
    % Update the line with the new sliced data
    for i = 1:length(curLineHandle)
        set(curLineHandle(i), 'xdata', get(curSlicer(i), 'XData'));
        set(curLineHandle(i), 'ydata', get(curSlicer(i), 'YData'));
        if isfield(geometry, 'Y')
            set(curLineHandle(i), 'zdata', get(curSlicer(i), 'ZData'));
            if isfield(geometry, 'Z')
                set(curLineHandle(i), 'cdata', get(curSlicer(i), 'CData'));
            end
        end
    end
    
    close(newfig);

    % Draw the plot!
    drawnow;
end

