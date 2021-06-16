function myslider

    x = 1:10;
    hplot = plot(x, 0*x);
    h = uicontrol('style', 'slider', 'units', 'pixel', 'position', [20 20 300 20]);
    addlistener(h, 'ContinuousValueChange', @(hObject, event) makeplot(hObject, event, x, hplot));
end

function makeplot(hObject, event, x, hplot)
    n = get(hObject, 'Value');
    set(hplot, 'ydata', x.^n);
    drawnow;
end
