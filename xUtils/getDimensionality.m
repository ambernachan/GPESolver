function [dimensions] = getDimensionality(geometry)

    dimensions = 1;
    if isfield(geometry, 'Y')
        dimensions = dimensions + 1;
        if isfield(geometry, 'Z')
            dimensions = dimensions + 1;
        end
    end
    
end