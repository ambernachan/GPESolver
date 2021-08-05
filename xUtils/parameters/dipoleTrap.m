function [trap] = dipoleTrap(parameters, X,Y,Z)
    
    if strcmp(parameters.atom, 'Na')
        Udp0 = getsimconst('dipoleTrap0_Na');
    elseif strcmp(parameters.atom, 'Rb')
        Udp0 = getsimconst('dipoleTrap0_Rb');
    elseif strcmp(parameters.atom, 'ferro')
        Udp0 = getsimconst('dipoleTrap0_Rb');
    else
        error('atom type is not recognized with this input of parameters.')
    end
    
    w0x = getsimconst('dipole_waist_x');
    w0y = getsimconst('dipole_waist_y');
    wx = beam_width(w0x, getsimconst('zRx'));
    wy = beam_width(w0y, getsimconst('zRy'));
    
%     trap = cell(3,3);
%     for n = 1:parameters.nComponents
%         for m = 1:parameters.nComponents
%             if n == m % diagonal elements
%                 trap{n,m} = Udp0 * w0x*w0y/(wx(Z).*wy(Z)) ...
%                     .* exp(-2*X.^2 ./ wx(Z).^2) .* exp(-2*Y.^2 ./ wy(Z).^2);
%             else
%                 trap{n,m} = 0;
%             end
%         end
%     end
    
    trap = Udp0 * w0x*w0y/(wx(Z).*wy(Z)) ...
        .* exp(-2*X.^2 ./ wx(Z).^2) .* exp(-2*Y.^2 ./ wy(Z).^2);
    
end