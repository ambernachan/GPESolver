
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
    end
    methods
        % Constructor
        function obj = PhiData(phi, geometry)

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

            if ~iscell(phi)
                    obj.phi = {phi};
            else
                obj.phi = phi;
            end
           
            obj.ncomponents = size(phi,2);
            
            % create obj.phisq, obj.phase
            obj.make_phisq_and_phase_from_phi();
            
            % normalizes phi using phisq, then redefines obj.phisq, obj.phase
            obj.normalize(); 

        end

        % function L2norm output: computation of the L2 norm in 1,2,3 dimensions,
            %   L2norm1d = sqrt(obj.dx) * sqrt(sum(obj.phisq));
            %   L2norm2d = sqrt(obj.dx*obj.dy) * sqrt(sum(sum(obj.phisq)));
            %   L2norm3d = sqrt(obj.dx*obj.dy*obj.dz) * sqrt(sum(sum(sum(obj.phisq))));
        function [L2norm] = L2norm(obj, phisqarray)
            phisum = cell(1, obj.ncomponents);
            prefactor = 1;
            axisNames = ['x', 'y', 'z'];

            for dim = 1:obj.dimensions
                field = ['d' axisNames(dim)]; % field = 'dx', 'dy', 'dz'
                prefactor = prefactor * (obj.(field));

                % creating the sum term (see comment below)
                if dim == 1
                    phisum = sum(phisqarray);
                else
                    phisum = sum(phisum);
                end
            end
            L2norm = sqrt(prefactor) * sqrt(phisum);
        end

    end

    methods (Access = private)

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


    end

end
