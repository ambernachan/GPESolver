classdef ExpPhis1D < ExpPhisBase
    properties

    end
    methods
        % Constructor
        function obj = ExpPhis1D(chi, gammas, r0, geometry)
            X = geometry.X;
            dx = geometry.dx;
            obj = obj@ExpPhisBase(chi, gammas, r0, X, X, X, dx, dx, dx);
            obj.ncomponents = 1;
        end
        
    end
end

