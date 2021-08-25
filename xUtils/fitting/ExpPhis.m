classdef ExpPhis < handle
    properties %(Access = private)
        %ncomponents % haven't yet implemented this, would have to change beta/chi/w?
        dr

        X1d

        X2d
        Y2d

        X3d
        Y3d
        Z3d

        gamma % struct with fields x,y,z
        beta
        chi
        r0 % struct with fields x,y,z
    end
    methods
        % Constructor
        function obj = ExpPhis(chi, gammas, r0, geometry)
            if nargin ~= 4
                error('Not enough, or too many, arguments (%d). Provide the following: chi, gammas, r0 and geometry.', nargin)
            end

            [obj.X1d, obj.X2d, obj.Y2d, obj.X3d, obj.Y3d, obj.Z3d, obj.dr] = ExpPhis.parseGeometry(geometry);
            obj.r0 = struct('x', r0(1), 'y', r0(2), 'z', r0(3));
            obj.chi = chi;
            obj.gamma = struct('x', gammas(1), 'y', gammas(2), 'z', gammas(3));
            obj.beta = 4*pi*obj.chi;
        end % end of constructor

        function [gaussapprox,w] = Gaussian1d(obj)
            w = obj.find_1d_w();
            prefactor = obj.gamma.x^(1/4) / (pi^(1/4) * sqrt(w));
            exponent = - obj.gamma.x * (obj.X1d - obj.r0.x).^2 / (2*w^2);
            gaussapprox = prefactor * exp(exponent);

            geom = struct('X', obj.X1d, 'dx', obj.dr.x);
            gaussapprox = PhiData(gaussapprox, geom, 'phi');
        end

        function [gaussapprox,w] = Gaussian2d(obj)
            w = obj.find_2d_w();
            prefactor = (obj.gamma.x*obj.gamma.y)^(1/4) / (sqrt(pi) * w);
            exponentx = - obj.gamma.x * (obj.X2d - obj.r0.x).^2 / (2*w^2);
            exponenty = - obj.gamma.y * (obj.Y2d - obj.r0.y).^2 / (2*w^2);
            gaussapprox = prefactor .* exp(exponentx) .* exp(exponenty);

            geom = struct('X', obj.X2d, 'Y', obj.Y2d, 'dx', obj.dr.x, 'dy', obj.dr.y);
            gaussapprox = PhiData(gaussapprox, geom, 'phi');
        end

        function [gaussapprox, w] = Gaussian3d(obj)
            w = obj.find_3d_w();
            prefactor = (obj.gamma.x*obj.gamma.y*obj.gamma.z)^(1/4) / ( pi^(3/4) * w^(3/2) );
            exponentx = - obj.gamma.x * (obj.X3d - obj.r0.x).^2 / (2*w^2);
            exponenty = - obj.gamma.y * (obj.Y3d - obj.r0.y).^2 / (2*w^2);
            exponentz = - obj.gamma.z * (obj.Z3d - obj.r0.z).^2 / (2*w^2);
            gaussapprox = prefactor .* exp(exponentx) .* exp(exponenty) .* exp(exponentz);

            geom = struct('X', obj.X3d, 'Y', obj.Y3d, 'Z', obj.Z3d, 'dx', obj.dr.x, 'dy', obj.dr.y, 'dz', obj.dr.z);
            gaussapprox = PhiData(gaussapprox, geom, 'phi');
        end

        function [TFapprox,redge] = TF1d(obj)
            V = 0.5 * obj.gamma.x^2 * obj.X1d.^2;
            TFapprox = (0.5 * (3*obj.beta*obj.gamma.x/2)^(2/3) - V) / obj.beta;

            %redge = (6*pi*obj.chi*obj.gamma.x )^(1/3); % find the cloud edge
            redge = (6*pi*obj.chi/obj.gamma.x^2)^(1/3); % find the cloud edge
            R = obj.X1d;

            geom = struct('X', obj.X1d, 'dx', obj.dr.x);

            TFapprox = ExpPhis.ApplyBoundaryAndConvertToPhiData(TFapprox, redge, R, geom, obj.gamma);
        end

        function [TFapprox,redge] = TF2d(obj)
            V = 0.5 * obj.gamma.x^2 * obj.X2d.^2 + 0.5 * obj.gamma.y^2 * obj.Y2d.^2;
            TFapprox = (sqrt(obj.beta*obj.gamma.x*obj.gamma.y / pi) - V) / obj.beta;

            %redge = (16*obj.chi * obj.gamma.x * obj.gamma.y)^(1/4); % find the cloud edge
            %R = sqrt( obj.X2d.^2 + obj.Y2d.^2 );
            
            redge = (16*obj.chi * obj.gamma.x * obj.gamma.y)^(1/4); % find the cloud edge
            R = sqrt( obj.gamma.x^2*obj.X2d.^2 + obj.gamma.y^2*obj.Y2d.^2 );

            geom = struct('X', obj.X2d, 'Y', obj.Y2d, 'dx', obj.dr.x, 'dy', obj.dr.y);

            TFapprox = ExpPhis.ApplyBoundaryAndConvertToPhiData(TFapprox, redge, R, geom, obj.gamma);
            redge = TFapprox.edge.dim2;
        end

        function [TFapprox,redge] = TF3d(obj)
            V = 0.5 * obj.gamma.x^2 * obj.X3d.^2 + 0.5 * obj.gamma.y^2 * obj.Y3d.^2 + 0.5 * obj.gamma.z^2 * obj.Z3d.^2;
            TFapprox = (0.5 * (15*obj.beta*obj.gamma.x*obj.gamma.y*obj.gamma.z / (4*pi))^(2/5) - V) / obj.beta;

            %redge = (15*obj.chi * obj.gamma.x * obj.gamma.y * obj.gamma.z)^(1/5); % find the cloud edge
            %R = sqrt( obj.X3d.^2 + obj.Y3d.^2 + obj.Z3d.^2);
            
            % find the cloud edge
            redge = (15*obj.chi)^(1/5);
            R = sqrt( obj.gamma.x^2 * obj.X3d.^2 + obj.gamma.y^2 * obj.Y3d.^2 + obj.gamma.z^2 * obj.Z3d.^2);

            geom = struct('X', obj.X3d, 'Y', obj.Y3d, 'Z', obj.Z3d, 'dx', obj.dr.x, 'dy', obj.dr.y, 'dz', obj.dr.z);

            TFapprox = ExpPhis.ApplyBoundaryAndConvertToPhiData(TFapprox, redge, R, geom, obj.gamma);
            redge = TFapprox.edge.dim3;
        end

        function [w] = getW(obj)
            w = struct();
            w.d1 = obj.find_1d_w();
            w.d2 = obj.find_2d_w();
            w.d3 = obj.find_3d_w();
        end
        
        function [x3d,y3d,z3d] = getRarrays(obj)
            x3d=obj.X3d;
            y3d=obj.Y3d;
            z3d=obj.Z3d;
        end

    end % end of methods

    % private (only accesible by this class) methods.
    methods(Access = private)
        function [w1d] = find_1d_w(obj)
            syms w positive
            eqn = w^4 - w*obj.beta/sqrt(2*pi*obj.gamma.x) - 1 == 0;
            sol = solve(eqn, w, 'Real', true, 'MaxDegree', 4);
            w1d = double(sol);
            w1d = w1d(w1d>0);
        end

        function [w2d] = find_2d_w(obj)
            w2d = ( 1 + (4*obj.chi)*sqrt(obj.gamma.x * obj.gamma.y)/(obj.gamma.x + obj.gamma.y) )^(1/4);
        end

        function [w3d] = find_3d_w(obj)
            syms w positive
            eqn = w^5 - w - 3*obj.beta*sqrt(obj.gamma.x*obj.gamma.y*obj.gamma.z) / ( (2*pi)^(3/2) * (obj.gamma.x+obj.gamma.y+obj.gamma.z) ) == 0;
            w3d = vpasolve(eqn, w, [0 Inf]);
            w3d = double(w3d);
            w3d = w3d(w3d>0);
        end
    end

    % Static (do not require an instance of the class because they do not
    % change or require any properties of the class) and private (only
    % accessible by this class) methods.
    methods (Static, Access = private)
        % Apply boundaries to function F , sets all values of F outside boundary R > redge or R < -redge to zero
        % Creates phidata from result, setting geom and redge
        function [BoundedApproxPhiData] = ApplyBoundaryAndConvertToPhiData(F, redge, R, geom, gammas)
            BoundedF = ExpPhis.ApplyBoundary(F, redge, R);
            BoundedApproxPhiData = PhiData(BoundedF, geom, 'phisq');
            % Store edge in PhiData struct
            BoundedApproxPhiData.setEdge(redge,gammas);
        end

        function [BoundedF] = ApplyBoundary(F, redge, R)
            indexes = or(R > redge, R < -redge);
            F(indexes) = 0;
            BoundedF = F;
        end

        function [X1d, X2d, Y2d, X3d, Y3d, Z3d, dr] = parseGeometry(geometry)
            if isfield(geometry, 'X') && isfield(geometry, 'Y') && isfield(geometry, 'Z') % 3D
                X = geometry.X;
                Y = geometry.Y;
                Z = geometry.Z;
                dx = geometry.dx;
                dy = geometry.dy;
                dz = geometry.dz;

                xSize = size(X, 2); % in 3D meshgrid 2nd dimension is X
                ySize = size(Y, 1); % in 3D meshgrid 1st dimension is Y
                zSize = size(Z, 3); % in 3D meshgrid 3th dimension is Z
                Xrange = linspace(min(min(min(X))), max(max(max(X))), xSize);
                Yrange = linspace(min(min(min(Y))), max(max(max(Y))), ySize);
                Zrange = linspace(min(min(min(Z))), max(max(max(Z))), zSize);
            elseif isfield(geometry, 'X') && isfield(geometry, 'Y') %2D
                X = geometry.X;
                Y = geometry.Y;
                dx = geometry.dx;
                dy = geometry.dy;

                xSize = size(X, 2); % in 2D meshgrid 2nd dimension is X
                ySize = size(Y, 1); % in 2D meshgrid 1st dimension is Y
                xMin = min(min(X));
                xMax = max(max(X));
                yMin = min(min(Y));
                yMax = max(max(Y));
                % For z we pick the minimum and maximum of the x and y range and the maximum size of those aswell
                zMin = min(xMin, yMin);
                zMax = max(xMax, yMax);

                % Take as delta z the minimal delta of x and y and calculate
                % the corresponding size
                dz = min(dx, dy);
                zSize = ((zMax - zMin) / dz) + 1;

                Xrange = linspace(xMin, xMax, xSize);
                Yrange = linspace(yMin, yMax, ySize);
                Zrange = linspace(zMin, zMax, zSize);
            elseif isfield(geometry, 'X') %1D
                Xrange = geometry.X;
                Yrange = Xrange;
                Zrange = Xrange;
                dx = geometry.dx;
                dy = dx;
                dz = dx;
            else
                error('Invalid geometry, provide geometry as a struct with X, X & Y or X & Y & Z members and corresponding d{x,y,z} value(s).')
            end
            X1d = Xrange;
            [X2d, Y2d] = meshgrid(Xrange, Yrange);
            [X3d, Y3d, Z3d] = meshgrid(Xrange, Yrange, Zrange);
%             X3d = 1;
%             Y3d = 1;
%             Z3d = 1;
            dr = struct('x', dx, 'y', dy, 'z', dz);
        end
    end

end
