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

    % a = 1/2*beta and b = gammax gammay gammaz
    tfEqn = ' real(sqrt( a^(3/5)*(15*b/(8*pi))^(2/5) - a*x.^2 )).^2 ';
    guess.a = 1/(2*info.params.betan);
    guess.b = info.params.gammas(1) * info.params.gammas(2) * info.params.gammas(3);
    
    f = fit(xdata',ydata',tfEqn,'Normalize','off', 'Robust' , 'bisquare', 'StartPoint',[guess.a, guess.b]);
    
    %Get x location of max value
    coeffs = coeffvalues(f);
    A = coeffs(1);
    B = coeffs(2);
    
    uncertainties = confint(f);
    db3 = abs(coeffs(2)-uncertainties(3));
    db4 = abs(coeffs(2)-uncertainties(3));
    uncertainties(3) = B - db3*2*coeffs(2); %lower B
    uncertainties(4) = B + db4*2*coeffs(2); %upper B

%     % %Increase resolution of x data (by 30)
%     xdataFine=(linspace(xdata(1),xdata(end),30))';  
    % Create high res Gauss function
    modelTF =  real(sqrt( A^(3/5)*(15*B/(8*pi))^(2/5) - A*X{1}.^2 )).^2;
    %Find max value and its x location
    ypeak = max(modelTF);
    Redge = (15*B / (8*pi*A))^(1/5);
    
    ydata = ydata / L2_norm1d(ydata, geometry);
    modelTF = modelTF / L2_norm1d(modelTF, geometry);
    plot(X{1}, ydata)
    hold on;
    plot(X{1}, modelTF)
end