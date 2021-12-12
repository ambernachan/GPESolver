function [num] = elementIndex(nthElement, elementNames, multiplicity, currentRunNumber)
    
    % Sanity check
    maxRunNumber = 1;
    for i = 1:numel(multiplicity)
        maxRunNumber = maxRunNumber*multiplicity(i);
    end
    if currentRunNumber > maxRunNumber
        warning('Your current run number exceeds the allowed maximum')
        return;
    end
    
    rightlength = 1;
    for k = 1:numel(elementNames)
        if k <= nthElement
            continue;
        else
            rightlength = rightlength*multiplicity(k);
        end
    end
    
    if mod(floor((currentRunNumber - 1) / rightlength)+1, multiplicity(nthElement)) == 0
        num = multiplicity(nthElement);
    else
        num = mod(floor((currentRunNumber - 1) / rightlength)+1, multiplicity(nthElement));
    end
    
    sprintf('Element index #%d\nnth element: #%d (out of %d)\n#cr = %d / %d', num,...
        nthElement, numel(multiplicity), currentRunNumber, maxRunNumber)
    
end