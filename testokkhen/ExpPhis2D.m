
classdef ExpPhis2D < ExpPhisBase
    properties
    end
    methods
        % Constructor
        function obj = ExpPhis2D(chi, gammas, r0, geometry)
            X = geometry.X;
            Y = geometry.Y;
            dx = geometrty.dx;
            dy = geometrty.dy;
            Z = X;
            dz = geometrty.dx;
            if max(X) < max(y)
               Z = Y;
               dz = geometrty.dy;
            end
            obj = obj@ExpPhisBase(chi, gammas, r0, X, Y, Z, dx, dy, dz);
            %obj.ncomponents = 2;
        end
        
    end
end