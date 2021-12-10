function [direction, dir] = getDirection(direction)

    if ischar(direction)
        if strcmp(direction, 'x')
            dir = 1;
        elseif strcmp(direction, 'y')
            dir = 2;
        elseif strcmp(direction, 'z')
            dir = 3;
        end
    elseif direction == 1
        direction = 'x';
        dir = 1;
    elseif direction == 2
        direction = 'y';
        dir = 2;
    elseif direction == 3
        direction = 'z';
        dir = 3;
    else
        error('Something is wrong in specifying the direction of plotting');
    end
    
end