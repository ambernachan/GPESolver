function [A, B, uncertainties, ypeak, f, X, modelTF, Redge] = tf_fit(geometry, ydata, info, axis)
    
    params = info.params;
    
    if nargin > 3
        xax = axis;
    else
        xax = 'x';
    end
    
    dims = numel(params.boxlimits);
    
    if strcmp(xax, 'x')
        nx = 1;
    elseif strcmp(xax, 'y')
        nx = 2;
    elseif strcmp(xax, 'z')
        nx = 3;
    else
        error('Direction not specified properly, try "x", "y", or "z".')
    end
    
    N = params.Ngridpts;
    xMax = params.boxlimits(nx);
    xMin = -xMax;
    
    if iscell(ydata)
        if numel(ydata) > 1
            p = 0;
            m.Ncomponents = numel(ydata);
            ydata = normalize_global(m, geometry, ydata);
            for n = 1:numel(ydata)
                p = p + abs(ydata{n}).^2;
                
            end
            ydata = p;
        else
            ydata = abs(ydata{1}).^2;
        end
    else
        ydata = abs(ydata).^2;
    end
    
    if dims == 1
        ydata = ydata / L2_norm1d(ydata, geometry)^2;
    elseif dims == 2
        ydata = ydata / L2_norm2d(ydata, geometry)^2;
    elseif dims == 3
        ydata = ydata / L2_norm3d(ydata, geometry)^2;
    end
    
    dx = 2*xMax / (N-1);
    X{1} = xMin:dx:xMax;
    X{2} = meshgrid(X{1}, X{1});
    X{3} = meshgrid(X{1}, X{1}, X{1});
    
    alpha = [{'X'}, {'Y'}, {'Z'}];
    r = 0;
    for d = 1:dims
        r = r + params.gammas(d)^2 * geometry.(alpha{d}).^2;
    end
    r = sqrt(r);
    
    if dims == 3
        if strcmp(xax, 'x')
            xdata = r((N-1)/2, :, (N-1)/2);
            ydata = ydata((N-1)/2, :, (N-1)/2);
        elseif strcmp(xax, 'y')
            xdata = r(:, (N-1)/2, (N-1)/2);
            ydata = ydata(:, (N-1)/2, (N-1)/2);
        elseif strcmp(xax, 'z')
            xdata = r(:, (N-1)/2, (N-1)/2);
            ydata = ydata(:, (N-1)/2, (N-1)/2);
            xdata = reshape(xdata, [1, N]);
            ydata = reshape(ydata, [1, N]);
        end
    elseif dims == 1
        xdata = r;
    elseif dims == 2
        if strcmp(xax, 'x')
            xdata = r((N-1)/2, :);
            ydata = ydata((N-1)/2, :);
        elseif strcmp(xax, 'y')
            xdata = r(:, (N-1)/2);
            ydata = ydata(:, (N-1)/2);
        end
    end

    % a = 1/2*beta and b = gammax gammay gammaz
    guess.a = 1/(2*info.params.betan);
    guess.a = 1/(2*25064);
    guess.b = info.params.gammas(1);
    guess.b = 1;
    if dims == 1 % the 0.25 is Geometry1D.dx
        tfEqn = ' real(sqrt( a^(1/3)*(3*b/4)^(2/3) - a*x.^2 )).^2 / sqrt(0.25)*sqrt(sum(abs(real(sqrt( a^(1/3)*(3*b/4)^(2/3) - a*x.^2 )).^2).^2)) ';
%         tfEqn = ' real(sqrt( a^(1/3)*(3*b/4)^(2/3) - a*x.^2 )).^2 ';
    elseif dims == 3
        guess.b = guess.b * info.params.gammas(2) * info.params.gammas(3);
        tfEqn = ' real(sqrt( a^(3/5)*(15*b/(8*pi))^(2/5) - a*x.^2 )).^2 ';
    else
        guess.b = guess.b * info.params.gammas(2);
        tfEqn = ' real(sqrt( sqrt(a)*sqrt(2*b/pi) - a*x.^2 )).^2 ';
%         error('dims==2 not supported yet.')
    end
    
    f = fit(xdata,ydata,tfEqn, 'Normalize','off', 'Robust' , 'bisquare', 'StartPoint', [guess.a, guess.b]);
%     f = fit(xdata,ydata,tfEqn, 'Normalize','off', 'Robust' , 'bisquare', 'StartPoint', [guess.a]);
    
    %Get x location of max value
    coeffs = coeffvalues(f);
    A = coeffs(1);
    B = coeffs(2);
%     B=1;
    uncertainties = confint(f);
    db3 = abs(coeffs(2)-uncertainties(3));
    db4 = abs(coeffs(2)-uncertainties(3));
    uncertainties(3) = B - db3*2*coeffs(2); %lower B
    uncertainties(4) = B + db4*2*coeffs(2); %upper B

%     % %Increase resolution of x data (by 30)
%     xdataFine=(linspace(xdata(1),xdata(end),30))';  
    % Create high res Gauss function
    % Not sure whether the Redge-s should include gammas...
    if dims == 1
        modelTF =  real(sqrt( A^(1/3)*(3*B/4)^(2/3) - A*X{1}.^2 )).^2 / (sqrt(0.25)*sqrt(sum(abs(real(sqrt( A^(1/3)*(3*B/4)^(2/3) - A*X{1}.^2 )).^2).^2))) ;
        Redge = (3*B/(4*A))^(1/3);
    elseif dims == 2
        real(sqrt( sqrt(A)*sqrt(2*B/pi) - A*X{1}.^2 )).^2;
        Redge = (2*B / (pi*A))^(1/4);
    elseif dims == 3
        modelTF =  real(sqrt( A^(3/5)*(15*B/(8*pi))^(2/5) - A*X{1}.^2 )).^2;
        Redge = (15*B / (8*pi*A))^(1/5);
    end
    % Find max value and its x location
    ypeak = max(modelTF);
    ydata = ydata / L2_norm1d(ydata, geometry);
    modelTF = modelTF / L2_norm1d(modelTF, geometry);
    plot(X{1}, ydata)
    hold on;
    plot(X{1}, modelTF)
end