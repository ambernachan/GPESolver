function compareGTF(geometry, phi, params)
    
    dims = 1;
    N{1} = geometry.Nx;
    xMax = max(geometry.X, [], 'all');
    xMin = min(geometry.X, [], 'all');
    if isfield(geometry, 'dy')
        dims = dims + 1;
        N{2} = geometry.Ny;        
        if isfield(geometry, 'dz')
            dims = dims + 1;
            N{3} = geometry.Nz;
        end
    end
    
    if iscell(phi)
        if numel(phi) > 1
            p = 0;
            m.Ncomponents = numel(phi);
            phi = normalize_global(m, geometry, phi);
            reqnorm = false;
            for n = 1:numel(phi)
                p = p + abs(phi{n});
                
            end
            phi = p;
        else
            reqnorm = true;
            phi = abs(phi{1});
        end
    else
        reqnorm = true;
        phi = abs(phi);
    end
    
    if reqnorm
        if dims == 1
            phi = phi / L2_norm1d(phi, geometry);
        elseif dims == 2
            phi = phi / L2_norm2d(phi, geometry);
        elseif dims == 3
            phi = phi / L2_norm3d(phi, geometry);
        end
    end
    
    X{1} = xMin:geometry.dx:xMax;
    X{2} = meshgrid(X{1}, X{1});
    X{3} = meshgrid(X{1}, X{1}, X{1});
    
    alpha = [{'X'}, {'Y'}, {'Z'}];
    r = 0;
    for d = 1:dims
        r = r + params.gammas(d)^2 * geometry.(alpha{d}).^2;
    end
    r = sqrt(r);
    
    %% Creating Thomas-Fermi dist
    
    if dims == 1
        TF = ((3*params.betan*params.gammas(1)/2)^(2/3) - r.^2) / (2*params.betan);
        redge = 1;
    elseif dims == 2
        TF = sqrt(params.gammas(1)*params.gammas(2)/(params.betan*pi)) - 0.5*r.^2/params.betan;
        C = (4*params.betan / pi)^(1/4) / params.gammas(2);
        S = (params.gammas(1)/params.gammas(2))^2 - 1;
        if C == 0
            if S >= 0
                redge = 0.5*sqrt(S*(-1+sqrt(3)));
            elseif S < 0
                redge = 0.5*sqrt(S*(-1-sqrt(3)));
            end
            return;
        end
        if S == 0
            redge = C;
        elseif S > 0 % meaning gx > gy and gy = 1
            redge = C;
            error('unknown redge')
        elseif S < 0 % meaning gx < gy and gx = 1
            redge = C;
            error('unknown redge')
        end
    elseif dims == 3
        TF = ((15*params.betan/(4*pi))^(2/5) - r.^2) / (2*params.betan);
        if params.gammas(1) == params.gammas(2) && params.gammas(1) == params.gammas(3)
            redge = (15*params.betan/(4*pi))^(1/5);
        else
            redge = (15*params.betan/(4*pi))^(1/5);
            error('unknown redge')
        end
    end
    
    indexes = or(r>redge,r<-redge);
    TF(indexes) = 0;
    TF = sqrt(TF);
    
    %% Creating Gaussian
    
    if dims == 1
        syms w positive
        eqn = w^4 - w*params.betan / sqrt(2*pi*params.gammas(1)) - 1 == 0;
        sol = solve(eqn, w, 'Real', true, 'MaxDegree', 4);
        W = double(sol);
        W = W(W>0);
    elseif dims == 2
        W = ( 1 + (4*params.chin)*sqrt(params.gammas(1) * params.gammas(2))/(params.gammas(1)+params.gammas(2)) )^(1/4);
    elseif dims == 3
        syms w positive
        eqn = w^5 - w - 3 * params.betan * sqrt( params.gammas(1)* ...
            params.gammas(2)*params.gammas(3) ) / ( (2*pi)^(3/2) * ...
            (params.gammas(1)+params.gammas(2)+params.gammas(3)) ) == 0;
        W = vpasolve(eqn, w, [0 Inf]);
        W = double(W);
        W = W(W>0);
    end
    
    GAUSS = 1;
    for d = 1:dims
        GAUSS = GAUSS * ( params.gammas(d)/(pi*W^2) )^(1/4);
    end
    GAUSS = GAUSS * exp(-r.^2 / (2*W^2));
    
    % Renormalizing everything
    if dims == 1
        phi = phi / L2_norm1d(phi, geometry);
        TF = TF / L2_norm1d(TF, geometry);
        GAUSS = GAUSS / L2_norm1d(GAUSS, geometry);
    elseif dims == 2
        phi = phi / L2_norm2d(phi, geometry);
        TF = TF / L2_norm2d(TF, geometry);
        GAUSS = GAUSS / L2_norm2d(GAUSS, geometry);
    elseif dims == 3
        phi = phi / L2_norm3d(phi, geometry);
        TF = TF / L2_norm3d(TF, geometry);
        GAUSS = GAUSS / L2_norm3d(GAUSS, geometry);
    end
    
    plot(X{1}, phi((N{2}-1)/2,:,(N{3}-1)/2));
    hold on;
    plot(X{1}, TF((N{2}-1)/2,:,(N{3}-1)/2));
    plot(X{1}, GAUSS((N{2}-1)/2,:,(N{3}-1)/2));
end