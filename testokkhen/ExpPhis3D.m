classdef ExpPhis3D < ExpPhisBase
    properties
    end
    methods
        % Constructor
        function obj = ExpPhis3D(chi, gammas, r0, geometry)
            X = geometry.X;
            Y = geometry.Y;
            Z = geometry.Z;
            dx = geometrty.dx;
            dy = geometrty.dy;
            dz = geometrty.dz;     
            obj = obj@ExpPhisBase(chi, gammas, r0, X, Y, Z, dx, dy, dz);
        end
        
    end
end