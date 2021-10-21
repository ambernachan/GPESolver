function compareGTF(geometry, phi, info, axis)
    
    params = info.params;
    
    if nargin > 3
        xax = axis;
    else
        xax = 'x';
    end
    
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
                p = p + abs(phi{n}).^2;
                
            end
            phi = p;
        else
            reqnorm = true;
            phi = abs(phi{1}).^2;
        end
    else
        reqnorm = true;
        phi = abs(phi).^2;
    end
    
    if reqnorm
        if dims == 1
            phi = phi / L2_norm1d(phi, geometry)^2;
        elseif dims == 2
            phi = phi / L2_norm2d(phi, geometry)^2;
        elseif dims == 3
            phi = phi / L2_norm3d(phi, geometry)^2;
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
        redge = (3*params.betan*params.gammas(1)/2)^(2/3);
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
    
    if dims == 3 % for |phi|^2
        TF = real(sqrt( (1/2*(15*params.betan*params.gammas(1)*params.gammas(2)*params.gammas(3)/(4*pi))^(2/5) - 0.5*r.^2)/params.betan )).^2;
    elseif dims == 1
        TF = real(sqrt( ((3*params.betan/2)^(2/3) - r.^2)/(2*params.betan) ));
    end
    
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
    
    GAUSS = 1; % for |phi|^2
    for d = 1:dims
        GAUSS = GAUSS * sqrt( params.gammas(d)/(pi*W^2) );
    end 
    GAUSS = GAUSS * exp(-r.^2 / (W^2));
    
%     % Renormalizing everything
%     if dims == 1
%         phi = phi / L2_norm1d(phi, geometry)^2;
%         TF = TF / L2_norm1d(TF, geometry)^2;
%         GAUSS = GAUSS / L2_norm1d(GAUSS, geometry)^2;
%     elseif dims == 2
%         phi = phi / L2_norm2d(phi, geometry)^2;
%         TF = TF / L2_norm2d(TF, geometry)^2;
%         GAUSS = GAUSS / L2_norm2d(GAUSS, geometry)^2;
%     elseif dims == 3
%         phi = phi / L2_norm3d(phi, geometry)^2;
%         TF = TF / L2_norm3d(TF, geometry)^2;
%         GAUSS = GAUSS / L2_norm3d(GAUSS, geometry)^2;
%     end

    if dims == 3
        TF = real(sqrt( (1/2*(15*params.betan*params.gammas(1)*params.gammas(2)*params.gammas(3)/(4*pi))^(2/5) - 0.5*X{1}.^2)/params.betan )).^2;
    end
    GAUSS = 1; % for |phi|^2
    for d = 1:dims
        GAUSS = GAUSS * sqrt( params.gammas(d)/(pi*W^2) );
    end 
    GAUSS = GAUSS * exp(-X{1}.^2 / (W^2));
    
    % Defining the data array for phi
    x = X{1};
    
    if dims == 1
        boxlim = params.boxlimits(1)
        phi_in = phi;
        tf_in = TF;
        gauss_in = GAUSS;
    elseif dims == 2
        error('2d plotting is not supported yet.')
    elseif dims == 3
        if strcmp(xax, 'x')
            boxlim = params.boxlimits(1);
    %         phi_in = phi((N{2}-1)/2, :, (N{3}-1)/2);
    %         tf_in = TF((N{2}-1)/2, :, (N{3}-1)/2);
    %         gauss_in = GAUSS((N{2}-1)/2, :, (N{3}-1)/2);
            if dims == 3
                phi_in = phi((N{2}-1)/2, :, (N{3}-1)/2);
                tf_in = TF;
                gauss_in = GAUSS;
            end
        elseif strcmp(xax, 'y')
            boxlim = params.boxlimits(2);
            phi_in = phi(:, (N{1}-1)/2, (N{3}-1)/2);
            tf_in = TF(:, (N{1}-1)/2, (N{3}-1)/2);
            gauss_in = GAUSS(:, (N{1}-1)/2, (N{3}-1)/2);
        elseif strcmp(xax, 'z')
            boxlim = params.boxlimits(3);
            phi_in = phi((N{2}-1)/2, (N{1}-1)/2, :);
            tf_in = TF((N{2}-1)/2, (N{1}-1)/2, :);
            gauss_in = GAUSS(:, (N{1}-1)/2, (N{3}-1)/2);
            phi_in = reshape(phi_in, [1, N{3}]);
            tf_in = reshape(tf_in, [1, N{3}]);
            gauss_in = reshape(gauss_in, [1, N{3}]);
        end
    end
    
    % renormalize again???
    phi_in = phi_in / L2_norm1d(phi_in, geometry);
    tf_in = tf_in / L2_norm1d(tf_in, geometry);
    gauss_in = gauss_in / L2_norm1d(gauss_in, geometry);
    
    % Creating figure
    evomarker = floor(length(x)/N{1});
    
    fig = figure();
    plot(x, phi_in, '-o', ...
        'MarkerIndices', 1:evomarker:length(phi_in), 'LineWidth', 1, ...
        'MarkerSize', 2, 'Color', [0 0 0.2]);
    hold on;
    plot(x, tf_in, '-', ...
        'LineWidth', 1, 'Color', [0.737 0.3216 0.5294]);
    plot(x, gauss_in, '-', ...
        'LineWidth', 1, 'Color', [0 0.6 0]);
    
    if params.chin > 1
        highesty = max(max(max(phi_in), max(tf_in)), max(gauss_in));
    else
        highesty = max(max(phi_in), max(gauss_in));
    end
    lowesty = min(min(max(phi_in), max(tf_in)), max(gauss_in));
    
    if highesty < 0
        maxlimy = highesty / 1.1;
    else % highesty > 0
        maxlimy = highesty * 1.1;
    end
    if lowesty < 0
        minlimy = lowesty * 1.1;
    else % lowesty > 0
        minlimy = min(min(abs(highesty)/100, 0), lowesty / 1.1);
    end
    
    ylim([minlimy maxlimy]);
    xlim([-boxlim, boxlim]);
    
    % Add axes labels and figure text
    xlabel([xax ' (a_{ho})']); 
    ylabel('|\phi|^2');
    title('Comparing |\phi| to gaussian and Thomas-Fermi distributions');
    
    lgd = legend('|\psi|^2', '|\psi_{TF}|^2', '|\psi_{gauss}|^2');
    
    savename = 'Compare G,TF to phi';
    
    %add datestring to figure
    annotation('textbox', [0, 0.05, 0, 0], 'string', sprintf('%s', info.creationTimeString))
    
    %add atom type text to figure
    atom_str = getAtomStr(params.atom);
    % add interaction parameters text to figure
    a = sprintf('%.2g', params.an); 
    ab = getphysconst('abohr'); nm = 10^(-9);
    a_ab = sprintf('%.2g', params.an/ab);
    a_nm = sprintf('%.2g', params.an/nm);
    n = sprintf('%1.4g', params.N); chi = sprintf('%1.4g', params.chin);
    param_str = {sprintf('%s', ['chi = ' chi]), sprintf('%s', ['at N = ' n ',']), ...
        sprintf('%s', ['   as = ' a_ab ' aB']), sprintf('%s', ['        = ' a_nm ' nm'])};
    
    % Add annotation about atom type to graph
    annotation('textbox', [0.15, 0.8, 0.1, 0.1], ...
        'string', sprintf('%s', atom_str), 'FitBoxToText', 'on')
    % Add annotation about interaction parameters to graph
    annotation('textbox', [0.15, 0.72, 0.1, 0.1], ...
        'string', param_str, 'FitBoxToText', 'on')
    
    % Save figure
    fig = gcf;
    info.save_figure(fig.Number, savename, '')
    info.save_figure(fig.Number, savename, '', info.fulldir, '.png')
    hold off
end