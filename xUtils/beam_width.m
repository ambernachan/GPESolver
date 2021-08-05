function [width] = beam_width(waist, zR)

    width = @(z) waist * sqrt(1+(z/zR).^2);
    
end