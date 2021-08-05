function [addedTraps] = addingPotentials(parameters, trap1, trap2)
    
    addedTraps = cell(3,3);
    for n = 1:parameters.nComponents
        for m = 1:parameters.nComponents
            addedTraps{n,m} = @(X,Y,Z) ( trap1 + trap2 );
        end
    end
    
end