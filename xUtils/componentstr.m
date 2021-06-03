function [component_string] = componentstr(nComponents)
    if nComponents == 3
        component_string = ['+', '0', '-'];
    elseif nComponents == 2
        component_string = ['+', '-'];
    elseif nComponents == 5
        component_string = ['+2','+1','0','-1','-2'];
    else
        component_string = '';
    end
end