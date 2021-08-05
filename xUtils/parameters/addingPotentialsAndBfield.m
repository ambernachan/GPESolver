function [addedTraps] = addingPotentialsAndBfield(parameters, trap, magnFieldPot)
    
    addedTraps = cell(3,3);
    for n = 1:parameters.nComponents
        for m = 1:parameters.nComponents
            mf = 0;
            if (n == m) && (n == 1)
                mf = 1;
            elseif (n == m) && (n == 3)
                mf = -1;
            end
            addedTraps{n,m} = @(X,Y,Z) ( @(X,Y,Z) trap + @(X,Y,Z) magnFieldPot(mf) );
        end
    end
    
end