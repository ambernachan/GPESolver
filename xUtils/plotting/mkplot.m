function mkplot(hObject, event, hplot, x, y)
    n = get(hObject, 'Value'); % y axis position (values in [0, 1])
    Ltime = length(x);
    if Ltime > 1
        time = round(n*(Ltime-1));
    else
        time = 1;
    end
%     set(hplot, 'ydata', y{time})
    set(hplot, 'ydata', x.^n);
    drawnow;
end