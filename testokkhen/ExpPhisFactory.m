classdef ExpPhisFactory
    methods(Static)
        function [ExpPhis] = Make(chi, gammas, r0, geometry)
            if nargin ~= 4
                error('Not enough, or too many, arguments, provide 4: chi, gammas, r0 and geometry.')
            end

            if isfield(geometry, 'X') && isfield(geometry, 'Y') && isfield(geometry, 'Z') % 3D
                ExpPhis = ExpPhis3D(chi, gammas, r0, geometry);
            elseif isfield(geometry, 'X') && isfield(geometry, 'Y') %2D
                ExpPhis = ExpPhis2D(chi, gammas, r0, geometry);
            elseif isfield(geometry, 'X') %1D
                ExpPhis = ExpPhis1D(chi, gammas, r0, geometry);
            else
                error('Invalid geometry, provide geometry as a struct with X, X & Y or X & Y & Z members and corresponding d{x,y,z} value.')
            end
        end
    end
end