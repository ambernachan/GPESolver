% FMGAUSSFIT Create/alter optimization OPTIONS structure.
%   [fitresult,..., rr] = fmgaussfit(xx,yy,zz) uses ZZ for the surface 
%   height. XX and YY are vectors or matrices defining the x and y 
%   components of a surface. If XX and YY are vectors, length(XX) = n and 
%   length(YY) = m, where [m,n] = size(Z). In this case, the vertices of the
%   surface faces are (XX(j), YY(i), ZZ(i,j)) triples. To create XX and YY 
%   matrices for arbitrary domains, use the meshgrid function. FMGAUSSFIT
%   uses the lsqcurvefit tool, and the OPTIMZATION TOOLBOX. The initial
%   guess for the gaussian is places at the maxima in the ZZ plane. The fit
%   is restricted to be in the span of XX and YY.
%   See:
%       http://en.wikipedia.org/wiki/Gaussian_function
%          
%   Examples:
%     To fit a 2D gaussian:
%       [fitresult, zfit, fiterr, zerr, resnorm, rr] =
%       fmgaussfit(xx,yy,zz);
%   See also SURF, OMPTMSET, LSQCURVEFIT, NLPARCI, NLPREDCI.

%   Copyright 2013, Nathan Orloff.

function [fitresult, zfit, fiterr, zerr, resnorm, rr] = fmgaussfit(xx,yy,zz)

    %% Condition the data
    [xData, yData, zData] = prepareSurfaceData( xx, yy, zz );
    xyData = {xData,yData};

    %% Set up the startpoint
    [amplitude, index] = max(zData);
    x0 = xData(index); % guess that the gauss center is at its maximum
    y0 = yData(index); % guess that the gauss center is at its maximum
    angl = 45; % angle in degrees.
    sigmay = 1;
    sigmax = 1;
    z0 = median(zData(:))-std(zData(:));
    xmax = max(xData)+2;
    ymax = max(yData)+2;
    xmin = min(xData)-2;
    ymin = min(yData)-2;

    %% Set up fittype and options.
    Lower = [0, 0, 0, 0, xmin, ymin, 0];
    Upper = [Inf, 180, Inf, Inf, xmax, ymax, Inf]; % angles greater than 90 are redundant
    StartPoint = [amplitude, angl, sigmax, sigmay, x0, y0, z0]; % [amplitude, sigmax, sigmaxy, sigmay, x0, y0, z0];

    tols = 1e-16;
    options = optimset('Algorithm','levenberg-marquardt',...
        'Display','off',...
        'MaxFunEvals',5e2,...
        'MaxIter',5e2,...
        'TolX',tols,...
        'TolFun',tols,...
        'TolCon',tols ,...
        'UseParallel','always');

    %% perform the fitting
    [fitresult,resnorm,residual] = ...
        lsqcurvefit(@gaussian2D,StartPoint,xyData,zData,Lower,Upper,options);
    [fiterr, zfit, zerr] = gaussian2Duncert(fitresult,residual,xyData);
    rr = rsquared(zData, zfit, zerr);
    zfit = reshape(zfit,size(zz));
    zerr = reshape(zerr,size(zz));

end

function rr = rsquared(z,zf,ze)
    % reduced chi-squared
    dz = z-zf;
    rr = 1./(numel(z)-8).*sum(dz.^2./ze.^2); % minus 8 because there are 7 fit parameters +1 (DOF)
end

% function z = gaussian2D(par,xy)
%     % compute 2D gaussian
%     z = par(7) + ...
%         par(1)*exp(-( ( (xy{1}-par(5)).*cosd(par(2))  + (xy{2}-par(6)).*sind(par(2)) ) ./par(3)).^2 ...
%                    -( ( -(xy{1}-par(5)).*sind(par(2)) + (xy{2}-par(6)).*cosd(par(2)) ) ./par(4)).^2 );
% end
function f = gaussian2D(parameters,XYZ,dim3)
    
    % % compute 2D gaussian V IDENTICAL TO COMMENTED ORIGINAL
    % f = t + ...
    %     amplitude*exp(-( ( (XYZ{1}-x0).*cosd(angl)  + (XYZ{2}-y0).*sind(angl) ) ./widthx).^2 ...
    %                -( ( -(XYZ{1}-x0).*sind(angl) + (XYZ{2}-y0).*cosd(angl) ) ./widthy).^2 );

    % widthx/y is defined as in exp(-X^2/width^2) so 2*sigma^2=width^2
    % t is a constant offset
    [amplitude, angl, widthx, widthy, x0, y0, t] = parameters(:);
    dim = 2;
    if dim3
        dim = 3;
        [amplitude, angl, widthx, widthy, widthz, x0, y0, z0, t] = parameters(:);
    end

    % compute 2D gaussian
    Rotation_matrix = makehgtform('zrotate', angl);
    R = Rotation_matrix(1:dim, 1:dim); % make it 2d or 3d
    if dim == 2
        Var = [1/widthx^2 0; 0 1/widthy^2]; % variance matrix (2d)
    else
        Var = [1/widthx^2 0; 0 1/widthy^2; 0 1/widthz^2]; % variance matrix (3d)
    Var = R*(Var)*R.';

    f = t + amplitude ...
         * exp( -( Var(1,1)*(XYZ{1}-x0).^2 ) ) ...
        .* exp( -( Var(2,2)*(XYZ{2}-y0).^2 ) ) ...
        .* exp( -( (Var(1,2)+Var(2,1)) .* (XYZ{1}-y0) .* (XYZ{2}-y0) ) );

    if dim == 3
        f = f .* exp( -( Var(3,3)*(XYZ{3}-z0).^2 ) );
    end

end

function [dpar,zf,dzf] = gaussian2Duncert(par,resid,xy)
    % get the confidence intervals
    J = guassian2DJacobian(par,xy);
    parci = nlparci(par,resid,'Jacobian',J);
    dpar = (diff(parci,[],2)./2)';
    [zf,dzf] = nlpredci(@gaussian2D,xy,par,resid,'Jacobian',J);
end

function J = guassian2DJacobian(parameters,xy)
    [amplitude, angl, widthx, widthy, x0, y0, t] = parameters(:);
    % compute the jacobian
    x = xy{1}; y = xy{2};
    J(:,1) = exp(- (cosd(angl).*(x - x0) + sind(angl).*(y - y0)).^2./widthx.^2 - (cosd(angl).*(y - y0) - sind(angl).*(x - x0)).^2./widthy.^2);
    J(:,2) = -amplitude.*exp(- (cosd(angl).*(x - x0) + sind(angl).*(y - y0)).^2./widthx.^2 - (cosd(angl).*(y - y0) - sind(angl).*(x - x0)).^2./widthy.^2).*((2.*(cosd(angl).*(x - x0) + sind(angl).*(y - y0)).*(cosd(angl).*(y - y0) - sind(angl).*(x - x0)))./widthx.^2 - (2.*(cosd(angl).*(x - x0) + sind(angl).*(y - y0)).*(cosd(angl).*(y - y0) - sind(angl).*(x - x0)))./widthy.^2);
    J(:,3) = (2.*amplitude.*exp(- (cosd(angl).*(x - x0) + sind(angl).*(y - y0)).^2./widthx.^2 - (cosd(angl).*(y - y0) - sind(angl).*(x - x0)).^2./widthy.^2).*(cosd(angl).*(x - x0) + sind(angl).*(y - y0)).^2)./widthx^3;
    J(:,4) = (2.*amplitude.*exp(- (cosd(angl).*(x - x0) + sind(angl).*(y - y0)).^2./widthx.^2 - (cosd(angl).*(y - y0) - sind(angl).*(x - x0)).^2./widthy.^2).*(cosd(angl).*(y - y0) - sind(angl).*(x - x0)).^2)./widthy^3;
    J(:,5) = amplitude.*exp(- (cosd(angl).*(x - x0) + sind(angl).*(y - y0)).^2./widthx.^2 - (cosd(angl).*(y - y0) - sind(angl).*(x - x0)).^2./widthy.^2).*((2.*cosd(angl).*(cosd(angl).*(x - x0) + sind(angl).*(y - y0)))./widthx.^2 - (2.*sind(angl).*(cosd(angl).*(y - y0) - sind(angl).*(x - x0)))./widthy.^2);
    J(:,6) = amplitude.*exp(- (cosd(angl).*(x - x0) + sind(angl).*(y - y0)).^2./widthx.^2 - (cosd(angl).*(y - y0) - sind(angl).*(x - x0)).^2./widthy.^2).*((2.*cosd(angl).*(cosd(angl).*(y - y0) - sind(angl).*(x - x0)))./widthy.^2 + (2.*sind(angl).*(cosd(angl).*(x - x0) + sind(angl).*(y - y0)))./widthx.^2);
    J(:,7) = ones(size(x));
end