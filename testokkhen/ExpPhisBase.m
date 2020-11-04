classdef ExpPhisBase < handle
    properties
    end
    properties(Access = private)
        %ncomponents
        dr

        TF
        Gaussian

        X1d
        
        X2d
        Y2d
        
        X3d
        Y3d
        Z3d

        gamma % struct with fields x,y,z
        beta
        chi
        r0
        w
    end
    methods
        % Constructor
        function obj = ExpPhisBase(chi, gammas, r0, X, Y, Z, dx, dy, dz)
            obj.chi = chi;
            obj.gamma = struct('x', gammas(1), 'y', gammas(2), 'z', gammas(3));
            obj.beta = 4*pi*obj.chi;
            obj.r0 = struct('x', r0(1), 'y', r0(2), 'z', r0(3));
            
            %obj.X = X;
            %obj.Y = Y;
            %obj.Z = Z;
            obj.dr.x = dx;
            obj.dr.y = dy;
            obj.dr.z = dz;
            
            obj.X1d = meshgrid(X);
            [obj.X2d, obj.Y2d] = meshgrid(X, Y);
            [obj.X3d, obj.Y3d, obj.Z3d] = meshgrid(X, Y, Z);
        end        
        
        function [gaussapprox] = Gaussian1d(obj)
            w = obj.find_1d_w();
            %obj.w.d1 = w;
            prefactor = obj.gamma.x^(1/4) / (pi^(1/4) * sqrt(w));
            exponent = - obj.gamma.x * (obj.X1d - obj.r0.x).^2 / (2*w^2);
            gaussapprox = prefactor * exp(exponent);

            geom = struct('X', obj.X1d, 'dx', obj.dr.x);
            gaussapprox = PhiData(gaussapprox, geom, 'phisq');

            %obj.Gaussian_1d = gaussapprox;
        end
                
        function [gaussapprox] = Gaussian2d(obj)
            w = obj.find_2d_w();
            %obj.w.d2 = w;
            prefactor = (obj.gamma.x*obj.gamma.y)^(1/4) / (sqrt(pi) * w);
            exponentx = - obj.gamma.x * (obj.X2d - obj.r0.x).^2 / (2*w^2);
            exponenty = - obj.gamma.y * (obj.Y2d - obj.r0.y).^2 / (2*w^2);
            gaussapprox = prefactor .* exp(exponentx) .* exp(exponenty);

            geom = struct('X', obj.X2d, 'Y', obj.Y2d, 'dx', obj.dr.x, 'dy', obj.dr.y);
            gaussapprox = PhiData(gaussapprox, geom, 'phisq');

            %obj.Gaussian_2d = gaussapprox;
        end
        
                
        function gaussian3D(obj)
            
        end
    end
    
    methods (Access = private)

        function [w1d] = find_1d_w(obj)
            syms w positive
            eqn = w^4 - w*obj.beta/sqrt(2*pi*obj.gamma.x) - 1 == 0;
            sol = solve(eqn, w, 'Real', true, 'MaxDegree', 4);
            w1d = double(sol)
        end

        function [w2d] = find_2d_w(obj)
            w2d = ( 1 + (4*obj.chi)*sqrt(obj.gamma.x * obj.gamma.y)/(obj.gamma.x + obj.gamma.y) )^(1/4);
        end

        function [w3d] = find_3d_w(obj)
            syms w positive
            eqn = w^5 - w - 3*obj.beta*sqrt(obj.gamma.x*obj.gamma.y*obj.gamma.z) / ( (2*pi)^(3/2) * (obj.gamma.x+obj.gamma.y+obj.gamma.z) ) == 0;
            w3d = vpasolve(eqn, w, [0 Inf])
            w3d = double(w3d)
        end

        % Apply boundaries to function F , sets all values of F outside boundary R > redge or R < -redge to zero
        % Creates phidata from result, setting geom and redge
        function [ApproxPhiData] = ApplyBoundaryAndConvertToPhiData(obj, F, redge, R, geom)
            indexes = find(or(R > redge, R < -redge));
            F(indexes) = 0;

            ApproxPhiData = PhiData(TFapprox, geom, 'phisq');
            ApproxPhiData.setEdge(redge); % set edge in PhiData struct
        end

    end
end

