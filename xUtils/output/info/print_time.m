function returnString = print_time(tElapsed)

% elapsed time format
if tElapsed < 60 % < 1 hour
    returnString = sprintf('%.4f seconds\n', tElapsed);

elseif tElapsed < (24*3600) % < 1 day
    hours = floor(tElapsed / 3600);
    minutes = floor( (tElapsed - hours*3600) / 60 );
    seconds = round(tElapsed - hours*3600 - minutes*60);
    returnString = sprintf('%02d:%02d:%02d (hh:mm:ss)\n', hours, minutes, seconds);

elseif tElapsed < (14*24*3600) % < 2 weeks
    days = floor(tElapsed / (3600*24));
    hours = floor( (tElapsed - days*(3600*24)) / 3600 );
    minutes = floor( (tElapsed - days*(3600*24) - hours*3600) / 60 );
    seconds = round(tElapsed - days*(3600*24) - hours*3600 - minutes*60);
    if days == 1
        returnString = sprintf('%d day & %02d:%02d:%02d (hh:mm:ss)\n', ...
            days, hours, minutes, seconds);
    else
        returnString = sprintf('%d days & %02d:%02d:%02d (hh:mm:ss)\n', ...
            days, hours, minutes, seconds);
    end
    
else % tElapsed > 2 weeks
    weeks = floor(tElapsed / (3600*24*7));
    days = floor( (tElapsed - weeks*(3600*24*7)) / (3600*24));
    hours = floor( (tElapsed - weeks*(3600*24*7) - days*(3600*24)) / 3600 );
    minutes = floor( (tElapsed - weeks*(3600*24*7) - days*(3600*24) - hours*3600) / 60 );
    seconds = round(tElapsed - weeks*(3600*24*7) - days*(3600*24) - hours*3600 - minutes*60);
    
    if weeks == 1
        returnString = sprintf('%d week, %dd, %dh, %dm, %ds\n', weeks, days, hours, minutes, seconds);
    else % weeks > 1
        returnString = sprintf('%d weeks, %dd, %dh, %dm, %ds\n', weeks, days, hours, minutes, seconds);
    end
    %returnString = strcat(returnString_a, '\n');
end

end