
classdef PhiData  < handle
    properties
        dimensions
        ncomponents

        phi
        phisq
        phase

        X
        Y
        Z
        dx
        dy
        dz

        edge
    end
    methods
        % Constructor
        function obj = PhiData(phi, geometry, phi_type)

            if ~iscell(phi)
                obj.phi = {phi};
            else
                obj.phi = phi;
            end

            obj.ncomponents = size(obj.phi,2);

            if nargin > 2
                if strcmp(phi_type, 'phisq')
                    obj.phisq = obj.phi;
                    obj.make_phi_and_phase_from_phisq();
                end
            end

            if nargin < 2 % geometry is not given
                obj.dimensions = ndims(obj.phi{1});
                geometry = struct();
                phi_type = input('Type either "phi"/"phisq" to define the data entered\n');
                if strcmp(phi_type, 'phisq')
                    obj.phisq = obj.phi; % def |phi|^2
                    obj.make_phi_and_phase_from_phisq(); % def phi, phase
                end
                
                box_size = input('What is the box size? ([-#x,#x]*[-#y,#y]*[-#z,#z], what is [#x #y #z]?)');
                if isempty(box_size)
                    box_size = [10 10 10];
                end
                
                for d = 1:obj.dimensions
                    if d == 1
                        xx = linspace(-box_size(d), box_size(d), size(obj.phi{1}, d));
                        geometry.dx = (max(xx)-min(xx)) / (size(obj.phi{1}, d) - 1);
                    end
                    if d == 2
                        yy = linspace(-box_size(d), box_size(d), size(obj.phi{1}, d));
                        geometry.dy = (max(yy)-min(yy)) / (size(obj.phi{1}, d) - 1);
                    end
                    if d == 3
                        zz = linspace(-box_size(d), box_size(d), size(obj.phi{1}, d));
                        geometry.dz = (max(zz)-min(zz)) / (size(obj.phi{1}, d) - 1);
                    end
                end

                if obj.dimensions == 1
                    geometry.X = xx;
                elseif obj.dimensions == 2
                    [geometry.X, geometry.Y] = meshgrid(xx,yy);
                elseif obj.dimensions == 3
                    [geometry.X, geometry.Y, geometry.Z] = meshgrid(xx,yy,zz);
                end
            end
                
            % Set number of dimensions, linspace X(,Y,Z) and grid spacing dx(,dy,dz)
            obj.dimensions = 1;
            obj.X = geometry.X;
            obj.dx = geometry.dx;
            if isfield(geometry, 'Y')
                obj.dimensions = obj.dimensions + 1;
                obj.Y = geometry.Y;
                obj.dy = geometry.dy;
                if isfield(geometry, 'Z')
                    obj.dimensions = obj.dimensions + 1;
                    obj.Z = geometry.Z;
                    obj.dz = geometry.dz;
                end
            end


            % create obj.phisq, obj.phase
            obj.make_phisq_and_phase_from_phi();
            
            % normalizes phi using phisq, then redefines obj.phisq, obj.phase
            obj.normalize(); 

        end % end of constructor

        function [geom] = return_geometry(obj)
            geom = struct('X',obj.X, 'dx',obj.dx);
            if obj.dimensions > 1
                geom.Y = obj.Y;
                geom.dy = obj.dy;
                if obj.dimensions > 2
                    geom.Z = obj.Z;
                    geom.dz = obj.dz;
                end
            end
            geom.ncomponents = obj.ncomponents;
            geom.dimensions = obj.dimensions;
        end
        
        function setEdge(obj, R, gammas)
            edge1d = R;
            edge2d = [R/gammas.x R/gammas.y];
            edge3d = [R/gammas.x R/gammas.y R/gammas.z];
            edge = struct();
            edge.dim1 = edge1d;
            edge.dim2 = edge2d;
            edge.dim3 = edge3d;
            obj.edge = edge;
        end

    end

    methods (Access = private)

        % function L2norm output: computation of the L2 norm in 1,2,3 dimensions,
        %   L2norm1d = sqrt(obj.dx) * sqrt(sum(obj.phisq));
        %   L2norm2d = sqrt(obj.dx*obj.dy) * sqrt(sum(sum(obj.phisq)));
        %   L2norm3d = sqrt(obj.dx*obj.dy*obj.dz) * sqrt(sum(sum(sum(obj.phisq))));
        function [L2norm] = L2norm(obj, phisqarray)
            % phisum = cell(1, obj.ncomponents); % this operation is per component n!
            phisqsum = [];
            prefactor = 1;
            axisNames = ['x', 'y', 'z'];

            for dim = 1:obj.dimensions
                field = ['d' axisNames(dim)]; % field = 'dx', 'dy', 'dz'
                prefactor = prefactor * (obj.(field));

                % creating the sum term (see comment below)
                if dim == 1
                    phisqsum = sum(phisqarray);
                else
                    phisqsum = sum(phisqsum);
                end
            end
            L2norm = sqrt(prefactor) * sqrt(phisqsum);
        end
        
        % normalizes phi using phisq, then redefines phi,phisq,phase
        function normalize(obj)
            L2norm = cell(1, obj.ncomponents);
            for n = 1:obj.ncomponents
                L2norm{n} = obj.L2norm(obj.phisq{n});
                obj.phi{n} = obj.phi{n} ./ L2norm{n};
            end
            obj.make_phisq_and_phase_from_phi();
        end
        
        % creates obj.phisq, obj.phase from obj.phi
        function make_phisq_and_phase_from_phi(obj)
            obj.phisq = cell(1, obj.ncomponents);
            obj.phase = cell(1, obj.ncomponents);
            for n = 1:obj.ncomponents
                obj.phisq{n} = abs(obj.phi{n}).^2;
                if isreal(obj.phi{n})
                    obj.phase{n} = [];
                else
                    obj.phase{n} = angle(obj.phi{n});
                end
            end
        end

        % creates obj.phisq, obj.phase from obj.phi
        function make_phi_and_phase_from_phisq(obj)
            obj.phi = cell(1, obj.ncomponents);
            obj.phase = cell(1, obj.ncomponents);
            for n = 1:obj.ncomponents
                obj.phi{n} = abs(sqrt(obj.phisq{n}));
                obj.phase{n} = [];
                %obj.phase{n} = phase(sqrt(obj.phisq{n})); % from |phi|^2 we can't reconstruct phase
            end
        end


    end

end
