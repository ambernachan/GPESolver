function [strS] = return_stringS(S)

    if mod(S,1) == 0 %checks whether S is an integer
        strS = sprintf('%d',S);
    elseif S >= 0.001
        strS = sprintf('%1.3f',S);
    elseif S < 0.001
        strS = sprintf('%.0e',S);
    end

end