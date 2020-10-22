
classdef ExpPhis  < handle
    properties
        ncomponents % haven't yet implemented this, would have to change beta/chi/w?
        dr

        TF_1d
        TF_2d
        TF_3d
        
        Gaussian_1d
        Gaussian_2d
        Gaussian_3d

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
        function obj = ExpPhis(chi, gammas, r0, geometry)

            if nargin == 0
                chi = input('Provide chi\n');

                gammas = input('Provide gammas [gx gy gz] or [] for gxyz = 1\n');
                if isempty(gammas)
                    gammas = [1, 1, 1];
                end
                
                r0 = input('Provide r0=[x0 y0 z0] for Gaussian f or [] for [000]/TF');
                if isempty(r0)
                    obj.r0 = struct('x', 0, 'y', 0, 'z', 0);
                else
                    x0 = r0(1);
                    if abs(r0(2)) > 0
                        y0 = r0(2);
                    else y0 = 0;
                    end
                    if abs(r0(3)) > 0
                        z0 = r0(3);
                    else z0 = 0;
                    end
                    obj.r0 = struct('x', x0, 'y', y0, 'z', z0);
                end
                
                gridxyz = input('Provide XYZ meshgrids [X Y Z]');
                dims = size(gridxyz,2) / size(gridxyz,1);
                if dims == 1
                    obj.X1d = gridxyz;
                    [obj.X2d, obj.Y2d] = meshgrid(obj.X1d, obj.X1d);
                    [obj.X3d, obj.Y3d, obj.Z3d] = meshgrid(obj.X1d, obj.X1d, obj.X1d);
                elseif dims == 2
                    obj.X2d = gridxyz(:, 1 : end/dims);
                    obj.Y2d = gridxyz(:, end/dims+1:end);
                    obj.X1d = obj.X2d(1,:);
                    y1d = obj.Y2d(1,:);
                    if max(y1d)>max(obj.X1d)
                        z1d = y1d;
                    else
                        z1d = obj.X1d;
                    end
                    [obj.X3d, obj.Y3d, obj.Z3d] = meshgrid(obj.X1d, y1d, z1d);
                elseif dims == 3
                    obj.X3d = gridxyz(:, 1 : end/dims, :);
                    obj.Y3d = gridxyz(:, end/dims+1:end*2/3, :);
                    obj.Z3d = gridxyz(:, end*2/3+1:end, :);
                    obj.X1d = obj.X3d(1,1:end,1);
                    y1d = obj.Y3d(1,1:end,1);
                    [obj.X2d, obj.Y2d] = meshgrid(obj.X1d, y1d);
                end

            end

            if nargin > 0
                obj.r0 = struct('x', r0(1), 'y', r0(2), 'z', r0(3));
            end
            
            if nargin < 4 % no geometry given
                box_size = input('What is the box size? ([-#x,#x]*[-#y,#y]*[-#z,#z], what is [#x #y #z]?)');
                Ngridpoints = input('What is the number of grid points per dimension? I.e. Nx=257');

                if size(box_size) == 1
                    box_size = [box_size(1), box_size(1), box_size(1)];
                elseif size(box_size) == 2
                    box_size = [box_size(1), box_size(2), max(box_size(1), box_size(2))];
                end
                
                xx = linspace(-box_size(d), box_size(d), Ngridpoints);
                yy = linspace(-box_size(d), box_size(d), Ngridpoints);
                zz = linspace(-box_size(d), box_size(d), Ngridpoints);
                
                obj.X1d = xx;
                [obj.X2d, obj.Y2d] = meshgrid(xx, yy);
                [obj.X3d, obj.Y3d, obj.Z3d] = meshgrid(xx, yy, zz);

                dx = (max(xx) - min(xx)) / (Ngridpoints - 1);
                dy = (max(yy) - min(yy)) / (Ngridpoints - 1);
                dz = (max(zz) - min(zz)) / (Ngridpoints - 1);

                obj.dr = struct('x', dx, 'y', dy, 'z', dz);

            elseif nargin == 4
                
                X = geometry.X;
                dx = geometry.dx;

                if isfield(geometry, 'Y') && ~isfield(geometry, 'Z') % 2d
                    Y = geometry.Y;
                    dy = geometry.dy;
                    
                    obj.X2d = X;
                    obj.Y2d = Y;
                    obj.X1d = X(1,:);
                    y1d = Y(1,:);
                    if max(X) > max(Y)
                        z1d = obj.X1d;
                        dz = dx;
                    else
                        z1d = y1d;
                        dz = dy;
                    end
                    [obj.X3d, obj.Y3d, obj.Z3d] = meshgrid(obj.X1d, y1d, z1d);
                    
                elseif isfield(geometry, 'Y') && isfield(geometry, 'Z') % 3d
                    Y = geometry.Y;
                    Z = geometry.Z;
                    dy = geometry.dy;
                    dz = geometry.dz;
                    
                    obj.X1d = X(1:length(X));
                    yy = Y(1:length(Y));
                    [obj.X2d, obj.Y2d] = meshgrid(obj.X1d, yy);
                    [obj.X3d, obj.Y3d, obj.Z3d] = meshgrid(X, Y, Z);
                    
                else % 1d
                    obj.X1d = X;
                    [obj.X2d, obj.Y2d] = meshgrid(X, X);
                    [obj.X3d, obj.Y3d, obj.Z3d] = meshgrid(X, X, X);
                    dy = geometry.dx;
                    dz = geometry.dx;
                end
                obj.dr = struct('x', dx, 'y', dy, 'z', dz);
            end

            obj.chi = chi;
            obj.gamma = struct('x', gammas(1), 'y', gammas(2), 'z', gammas(3));
            obj.beta = 4*pi*obj.chi;
            
            % w's only changed when function is called
            obj.w = struct('d1', 1, 'd2', 1, 'd3', 1);
            
        end % end of constructor

        function [gaussapprox] = Gaussian1d(obj)
            
            w = obj.find_1d_w();
            obj.w.d1 = w;
            prefactor = obj.gamma.x^(1/4) / (pi^(1/4) * sqrt(w));
            exponent = - obj.gamma.x * (obj.X1d - obj.r0.x).^2 / (2*w^2);
            gaussapprox = prefactor * exp(exponent);

            geom = struct('X', obj.X1d, 'dx', obj.dr.x);
            gaussapprox = PhiData(gaussapprox, geom, 'phisq');

            obj.Gaussian_1d = gaussapprox;
        end

        function [gaussapprox] = Gaussian2d(obj)
            
            w = obj.find_2d_w();
            obj.w.d2 = w;
            prefactor = (obj.gamma.x*obj.gamma.y)^(1/4) / (sqrt(pi) * w);
            exponentx = - obj.gamma.x * (obj.X2d - obj.r0.x).^2 / (2*w^2);
            exponenty = - obj.gamma.y * (obj.Y2d - obj.r0.y).^2 / (2*w^2);
            gaussapprox = prefactor .* exp(exponentx) .* exp(exponenty);

            geom = struct('X', obj.X2d, 'Y', obj.Y2d, 'dx', obj.dr.x, 'dy', obj.dr.y);
            gaussapprox = PhiData(gaussapprox, geom, 'phisq');

            obj.Gaussian_2d = gaussapprox;
        end

        function [gaussapprox] = Gaussian3d(obj)
            
            w = obj.find_3d_w();
            obj.w.d3 = w;
            prefactor = (obj.gamma.x*obj.gamma.y*obj.gamma.z)^(1/4) / ( pi^(3/4) * w^(3/2) );
            exponentx = - obj.gamma.x * (obj.X3d - obj.r0.x).^2 / (2*w^2);
            exponenty = - obj.gamma.y * (obj.Y3d - obj.r0.y).^2 / (2*w^2);
            exponentz = - obj.gamma.z * (obj.Z3d - obj.r0.z).^2 / (2*w^2);
            gaussapprox = prefactor .* exp(exponentx) .* exp(exponenty) .* exp(exponentz);

            geom = struct('X', obj.X3d, 'Y', obj.Y3d, 'Z', obj.Z3d, 'dx', obj.dr.x, 'dy', obj.dr.y, 'dz', obj.dr.z);
            gaussapprox = PhiData(gaussapprox, geom, 'phisq');

            obj.Gaussian_3d = gaussapprox;
        end

        function [TFapprox] = TF1d(obj)

            V = 0.5 * obj.gamma.x^2 * obj.X1d.^2;
            TFapprox = (0.5 * (3*obj.beta*obj.gamma.x/2)^(2/3) - V) / obj.beta;
            
            redge = (3*obj.beta / (2*obj.gamma.x) )^(1/3); % find the cloud edge
            sprintf('x_{edge} = %f',redge)

            indx = find(or(obj.X1d>redge,obj.X1d<-redge));
            TFapprox(indx) = 0; % remove data beyond the cloud edge

            geom = struct('X', obj.X1d, 'dx', obj.dr.x);
            TFapprox = PhiData(TFapprox, geom, 'phisq'); % construct PhiData

            TFapprox.redge(redge); % set Redge in PhiData struct

            obj.TF_1d = TFapprox;
        end

        function [TFapprox] = TF2d(obj)

            V = 0.5 * obj.gamma.x^2 * obj.X2d.^2 + 0.5 * obj.gamma.y^2 * obj.Y2d.^2;
            TFapprox = (sqrt(obj.beta*obj.gamma.x*obj.gamma.y / pi) - V) / obj.beta;

            geom = struct('X', obj.X2d, 'Y', obj.Y2d, 'dx', obj.dr.x, 'dy', obj.dr.y);
            TFapprox = PhiData(TFapprox, geom, 'phisq');
            
            obj.TF_2d = TFapprox;
        end

        function [TFapprox] = TF3d(obj)

            V = 0.5 * obj.gamma.x^2 * obj.X3d.^2 + 0.5 * obj.gamma.y^2 * obj.Y3d.^2 + 0.5 * obj.gamma.z^2 * obj.Z3d.^2;
            TFapprox = (0.5 * (15*obj.beta*obj.gamma.x*obj.gamma.y*obj.gamma.z / (4*pi))^(2/5) - V) / obj.beta;

            geom = struct('X', obj.X3d, 'Y', obj.Y3d, 'Z', obj.Z3d, 'dx', obj.dr.x, 'dy', obj.dr.y, 'dz', obj.dr.z);
            TFapprox = PhiData(TFapprox, geom, 'phisq');
            
            obj.TF_3d = TFapprox;
        end

    end % end of methods

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

    end

end
