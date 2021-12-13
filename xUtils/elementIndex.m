% elementIndex
%    num = elementIndex(nthElement, elementNames, multiplicity, 
%          currentRunNumber)
% Given a list of element names and the multiplicity of those elements
% (i.e. ["A1", "A2", "A3"] and [4, 5, 6], respectively), we can calculate a
% maxRunNumber (e.g. 4*5*6=120). Given a currentRunNumber we have a measure
% of which element in this (120-element-) series we are 'current'ly
% running. Given the nth element we know which element this is ("A1", "A2",
% or "A3"). The output of this function is 'num', which is the value of the
% current 'element' from the list of elements that currentRun is at. (E.g.
% in the as-above; if the nthElement = 3 and the currentRunNumber = x, we
% can find the how-manieth time we're running the simulation with element
% = 3, which can be 1-6 in this case).

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
    
%     rightlength is the value of length for the elements before nthElement
%     (E.g. if nthElement = 3 and multiplicity = [4, 5, 6] then rightlength
%     = 4*5 = 20).
    rightlength = 1;
    for k = 1:numel(elementNames)
        if k <= nthElement
            continue;
        else
            rightlength = rightlength*multiplicity(k);
        end
    end
    
%     if we're at the last value in the nth element, num =
%     multiplicity(nthElement), (e.g. if multiplicity = [4, 5, 6] and
%     nthElement = 3 and we're at the last value of element 3, num = 6).
    if mod(floor((currentRunNumber - 1) / rightlength)+1, multiplicity(nthElement)) == 0
        num = multiplicity(nthElement);
%     else if we're not at the last value; num = modulo[ floor{currRunNum-1
%     / rightlength}+1, multiplicity(nthElement) ], or: 
%     if multiplicity = [4, 5, 6] and nthElement = 3 and we're at the 2nd 
%     value of element 3, num = 2.
    else
        num = mod(floor((currentRunNumber - 1) / rightlength)+1, multiplicity(nthElement));
    end
    
    % Print current element info
    sprintf('Element index #%d\nnth element: #%d (out of %d)\n#cr = %d / %d', num,...
        nthElement, numel(multiplicity), currentRunNumber, maxRunNumber)
    
end