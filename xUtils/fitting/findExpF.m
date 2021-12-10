% gives the global, normalized <Fx>, <Fy>, <Fz> results in cell array
function [expectedF] = findExpF(Method, Geometry3D, Phi)
    
    expectedF = cell(1,3);
    
    % Gives normalized distribution of the Fx, Fy, Fz expectation values
    Fs = findExpFdistr(Method, Geometry3D, Phi);
    
    % Calculate the global <Fx>, <Fy>, <Fz>
    expectedF{1} = sum(sum(sum(Fs{1})));
    expectedF{2} = sum(sum(sum(Fs{2})));
    expectedF{3} = sum(sum(sum(Fs{3})));
    
end