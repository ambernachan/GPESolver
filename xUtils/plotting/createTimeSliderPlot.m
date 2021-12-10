function [sliderHandle] = createTimeSliderPlot(dataLength, figNumber)
%CREATETIMESLIDERPLOT Summary of this function goes here
%   Detailed explanation goes here

% Default: [0.3542 0.5167 0.2917 0.3889]
    if exist('figNumber', 'var')
        fig = figure(figNumber);
    else
        fig = figure();
    end
    fig.Units = 'normalized';
    fig.Position = [0.3542 0.5167 0.2917 0.3989];
    
    % Create the slider control
    dataRange = dataLength - 1;
    sliderHandle = uicontrol('style','slider',...
                'units','normalized',...
                'position',[0.937,0.102,0.032,0.824],...
                'min',1,'max',dataLength,'val',1,...
                'sliderstep',[1/dataRange 10/dataRange]);
              
end

