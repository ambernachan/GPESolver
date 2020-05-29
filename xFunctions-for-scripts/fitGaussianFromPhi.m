function fitGaussianFromPhi(phi, Geometry2D, varargin)

xbox = Geometry2D.X(end);
ybox = Geometry2D.Y(end);

%% Definition of the input options struct, with default values as given

% default struct: s = struct('fitmethod', 'both', 'centered', true, 'plotall', true,
% 'printall', true, 'showsource', true, 'dim3', false);

if nargin == 2
    s = OptionFactory();
elseif nargin == 3
    s = varargin{1};
elseif nargin > 3
    error('Too many input arguments')
else
    error('Something went wrong, check your input arguments')
end



%% Find array sizes

n = size(phi,1);
m = size(phi,2);

%% Find peak location and define corresponding data arrays 

[rowindex, colindex] = findMaxRowAndCol(phi);

if s.printall
    fprintf('Using data from x = %d/%d, y = %d/%d\n\n', colindex, n, rowindex, m);
end

xDATA = phi(rowindex, :);
yDATA = phi(:, colindex);

xa = linspace(-ybox,ybox,m);
ya = xDATA;
xb = linspace(-xbox,xbox,n);
yb = yDATA;

%% Derived quantities

if strcmp(s.fitmethod, 'both') || strcmp(s.fitmethod, 'manual')
    if ~s.centered
        [~, UexpectX] = max(xDATA);
        [~, UexpectY] = max(yDATA);
    else
        UexpectX = max(xa)/2;
        UexpectY = max(xb)/2;
    end
    WexpectX = max(xa)/8;
    WexpectY = max(xb)/8;
end

%% Generate (empty) output lists

widths1 = zeros(1, 6); % for gauss1 matlab std fit
widths = zeros(1, 3); % for gaussianFit output

%% View input image

if s.showsource
    fig=Figure_Var2d();
    fig.label = 1;
    fig.title = 'Input function \phi(x,y)';
    draw_function_2d(phi,Geometry2D, fig)
    
    %draw_phisq_2d(1, {phi}, Method, Geometry2D)
end

%%
%% FIT

%% gauss1 Matlab fit + saving parameters
% y = b*exp(-(x-x0)^2/w^2) fits coeffs b, x0, w

if strcmp(s.fitmethod, 'both') || strcmp(s.fitmethod, 'gauss1')
    
    fx = fit(xa.',ya.','gauss1');
    coeffx = coeffvalues(fx);
    fy = fit(xb.',yb,'gauss1');
    coeffy = coeffvalues(fy);

    gauss1CoeffsX = GaussParameterStruct(coeffx(1), coeffx(2), coeffx(3), confint(fx), 0);
    gauss1CoeffsY = GaussParameterStruct(coeffy(1), coeffy(2), coeffy(3), confint(fy), 0);
    
    % function coeffs = GaussParameterStruct(B, X0, W, varargin)
    %{
    arg4; varargin(1) = parameter A
    arg5; varargin(2) = confint(fitresult) for coefficients B, X0, W
    arg6; varargin(3) = confint for coefficient A
    %}
    % y = a + b*exp(-(x-x0)^2/w^2) fits coeffs A,B,U,W (=a, b, x0, w)
    
    % Determine average (x,y) widths/sigmas for 'gauss1'
    gauss1CoeffsAvg = avgGaussParameters(gauss1CoeffsX, gauss1CoeffsY);
    
end

%% Fitting using the script gaussianFit  + saving parameters
% y = a + b*exp(-(x-x0)^2/w^2) fits coeffs A,B,U,W (=a, b, x0, w)

if strcmp(s.fitmethod, 'both') || strcmp(s.fitmethod, 'manual')
    
    [Ax,Bx,Ux,Wx,confintX, xpeakx,ypeakx, ffx] = gaussianFit(xa.',ya.', UexpectX, WexpectX);
    [Ay,By,Uy,Wy,confintY, xpeaky,ypeaky, ffy] = gaussianFit(xb.',yb, UexpectY, WexpectY);

    gaussCoeffsX = GaussParameterStruct(Bx, Ux, Wx, confintX(3:8), Ax, confintX(1:2));
    gaussCoeffsY = GaussParameterStruct(By, Uy, Wy, confintY(3:8), Ay, confintY(1:2));
    
    % determine average (x,y) widths/sigmas for 'gaussianFit'
    gaussCoeffsAvg = avgGaussParameters(gaussCoeffsX, gaussCoeffsY);
    
end

%% Print outputs

if s.printall
    
    
    % For 'gauss1' default Matlab fit
    if strcmp(s.fitmethod, 'both') || strcmp(s.fitmethod, 'gauss1')
        
        % print fit results
        fprintf('gauss1 fit : y = b*exp(-(x-x0)^2/w^2), where w^2 = 2*sigma^2\n')
        fprintf('\t parameters (x) : a = X, b = %.4f+-%.4f, x0 = %.4f+-%.4f,\n\t\t\t\t w = %.4f+-%.4f and sigma = %.4f+-%.4f \n', ...
            gauss1CoeffsX.B, gauss1CoeffsX.uncB, gauss1CoeffsX.X0, gauss1CoeffsX.uncX0, ...
            gauss1CoeffsX.W, gauss1CoeffsX.uncW, gauss1CoeffsX.Sigma,gauss1CoeffsX.uncSigma)
        fprintf('\t parameters (y) : a = X, b = %.4f+-%.4f, x0 = %.4f+-%.4f,\n\t\t\t\t w = %.4f+-%.4f and sigma = %.4f+-%.4f \n', ...
            gauss1CoeffsY.B, gauss1CoeffsY.uncB, gauss1CoeffsY.X0, gauss1CoeffsY.uncX0, ...
            gauss1CoeffsY.W, gauss1CoeffsY.uncW, gauss1CoeffsY.Sigma,gauss1CoeffsY.uncSigma)
        
        % print widths
        fprintf('\t widths : \t wx = %.4f +- %.4f a_0, \twy = %.4f +- %.4f a_0. \n', ...
            gauss1CoeffsX.W, gauss1CoeffsX.uncW, gauss1CoeffsY.W, gauss1CoeffsY.uncW);
        fprintf('\t\t\t\t w ~ sqrt(wx wy) = %.3f +- %.3f a_0.\n', gauss1CoeffsAvg.W, gauss1CoeffsAvg.uncW);
        fprintf('\t std devs : sigmax = %.4f +- %.4f a_0, \t sigmay = %.4f +- %.4f a_0. \n', ...
            gauss1CoeffsX.Sigma, gauss1CoeffsX.uncSigma, gauss1CoeffsY.Sigma, gauss1CoeffsY.uncSigma);
        fprintf('\t\t\t\t sigma ~ sqrt(sigmax sigmay) = %.3f +- %.3f a_0.\n', gauss1CoeffsAvg.Sigma, gauss1CoeffsAvg.uncSigma);
    end
    
    % For 'gaussianFit' results
    if strcmp(s.fitmethod, 'both') || strcmp(s.fitmethod, 'manual')
        
        % print fit results
        fprintf('Gauss fit : y = a + b*exp(-(x-x0)^2/w^2), where w^2 = 2*sigma^2\n');
        fprintf('\t parameters (x) : a = %.4f+-%.4f, b = %.4f+-%.4f, x0 = %.4f+-%.4f,\n\t\t\t\t w = %.4f+-%.4f and sigma = %.4f+-%.4f \n', ...
            Ax, gaussCoeffsX.uncA, Bx, gaussCoeffsX.uncB, Ux, gaussCoeffsX.uncX0, ...
            Wx, gaussCoeffsX.uncW, gaussCoeffsX.Sigma,gaussCoeffsX.uncSigma);
        fprintf('\t parameters (y) : a = %.4f+-%.4f, b = %.4f+-%.4f, x0 = %.4f+-%.4f,\n\t\t\t\t w = %.4f+-%.4f and sigma = %.4f+-%.4f \n', ...
            Ay, gaussCoeffsY.uncA, By, gaussCoeffsY.uncB, Uy, gaussCoeffsY.uncX0, ...
            Wy, gaussCoeffsY.uncW, gaussCoeffsY.Sigma,gaussCoeffsY.uncSigma);
        
        % print widths
        fprintf('\t widths : \t wx = %.4f +- %.4f a_0, \t wy = %.4f +- %.4f a_0. \n', ...
            gaussCoeffsX.W, gaussCoeffsX.uncW, gaussCoeffsY.W, gaussCoeffsY.uncW)
        fprintf('\t\t\t\t w ~ sqrt(wx wy) = %.3f +- %.3f a_0.\n', gaussCoeffsAvg.W, gaussCoeffsAvg.uncW)
        fprintf('\t std devs : sigmax = %.4f +- %.4f a_0, \t sigmay = %.4f +- %.4f a_0. \n', ...
            gaussCoeffsX.Sigma, gaussCoeffsX.uncSigma, gaussCoeffsY.Sigma, gaussCoeffsY.uncSigma)
        fprintf('\t\t\t\t sigma ~ sqrt(sigmax sigmay) = %.3f +- %.3f a_0.\n', gaussCoeffsAvg.Sigma, gaussCoeffsAvg.uncSigma)
                
    end

end

%% X AND Y PLOT OF FITTED DATA

if s.plotall

    % to set figure size
    x0 = 0.1;
    y0 = 0.2;
    width = 0.8;
    height = 0.6;
    
    xlimit = [-xbox xbox];
    
    if strcmp(s.fitmethod, 'manual') || strcmp(s.fitmethod, 'both')
        xdataFine = (linspace(xa(1),xa(end),length(xa)))';
        fitGausx =  Ax + Bx * exp(- ((xdataFine-Ux).^2 ./ (Wx^2)) );
        xdataFine = (linspace(xa(1),xa(end),length(xa)))';

        ydataFine = (linspace(xb(1),xb(end),length(xb)))';  
        fitGausy =  Ay + By * exp(- ((ydataFine-Uy).^2 ./ (Wy^2)) );
        ydataFine = (linspace(xb(1),xb(end),length(xb)))';
        
        txtX = ['\leftarrow \sigma_x \approx ', sprintf('%.3f',gaussCoeffsX.Sigma), ...
            '+- ', sprintf('%.3f',gaussCoeffsX.uncSigma), ' [a_0]'];
        txtY = ['\leftarrow \sigma_y \approx ', sprintf('%.3f',gaussCoeffsY.Sigma), ...
            '+- ', sprintf('%.3f',gaussCoeffsY.uncSigma), ' [a_0]'];
        
        % limits
        ylimit = [0 max(max(fitGausx),max(ya))*1.1];
    end

    if strcmp(s.fitmethod, 'both') || strcmp(s.fitmethod, 'gauss1')
        txtX1 = ['\leftarrow \sigma_x \approx ', sprintf('%.3f',gauss1CoeffsX.Sigma), ...
            '+- ', sprintf('%.3f',gauss1CoeffsX.uncSigma), ' [a_0]'];
        txtY1 = ['\leftarrow \sigma_y \approx ', sprintf('%.3f',gauss1CoeffsY.Sigma), ...
            '+- ', sprintf('%.3f',gauss1CoeffsY.uncSigma), ' [a_0]'];
        
        % limits
        maxYdim1 = max(max(max(ya),max(yb)), max(fx.a1, fy.a1));
        ylimit1 = [0 maxYdim1*1.1];
    end
    
    if strcmp(s.fitmethod, 'both')
        ylimit = max(ylimit, ylimit1);
    elseif strcmp(s.fitmethod, 'gauss1')
        ylimit = ylimit1;
    end
    
    if strcmp(s.fitmethod, 'manual')
        figure(k)
        set(gcf,'units','normalized','position',[x0,y0,width,height])
        title('\phi');

        % X
        subplot(121)
        p11 = plot(xdataFine, fitGausx,'DisplayName','Fitted Data'); % fit data
        hold on;
        p12 = plot(xa, ya, 'DisplayName', 'Data'); % original data
        xlabel('x array [a_0]');
        ylabel('\phi(x)');
        title('fit (x)');
        text(0.75, 0.5, txtX, 'FontSize', 14, 'Units', 'normalized')
        hold off;
        xlim(xlimit);
        ylim(ylimit);
        legend([p11, p12])

        % Y
        subplot(122)
        p21 = plot(ydataFine, fitGausy); % fit data
        hold on;
        p22 = plot(xb, yb); % original data

        xlabel('y array [a_0]');
        ylabel('\phi(y)');
        title('fit (y)');
        text(0.75, 0.5, txtY, 'FontSize', 14, 'Units', 'normalized')
        hold off;
        xlim(xlimit);
        ylim(ylimit);
        legend([p21,p22]);
        
    end

    if strcmp(s.fitmethod, 'gauss1')
        figure(k)
        set(gcf,'units','normalized','position',[x0,y0,width,height])

        subplot(121)
        plot(fx, xa, ya) % data + fit

        % modify labels for tick marks
        xticks = get(gca,'xtick');
        newlabels = arrayfun(@(x) sprintf('%.1f', x), xticks, 'un', 0);
        set(gca,'xticklabel',newlabels);

        xlabel('[a_0]');
        ylabel('\phi_1(x)');
        title('gauss1 fit (x)');
        text(0.75, 0.5, txtX1, 'FontSize', 14, 'Units', 'normalized')
        
        xlim(xlimit);
        ylim(ylimit);
        
        subplot(122)
        plot(fy, xb, yb) % data + fit

        % modify labels for tick marks
        xticks = get(gca,'xtick');
        newlabels = arrayfun(@(x) sprintf('%.1f', x), xticks, 'un', 0);
        set(gca,'xticklabel',newlabels);

        xlabel('[a_0]');
        ylabel('\phi_1(y)');
        title('gauss1 fit (y)');
        text(0.75, 0.5, txtY1, 'FontSize', 14, 'Units', 'normalized')
        
        xlim(xlimit);
        ylim(ylimit);
        
    end

    if strcmp(s.fitmethod, 'both')
        figure()
        set(gcf,'units','normalized','position',[x0,y0,width,height])

        subplot(141)
        plot(fx, xa, ya) % data + fit
        
        % modify labels for tick marks
        xticks = get(gca,'xtick') ;
        newlabels = arrayfun(@(x) sprintf('%.1f', x), xticks, 'uni', 0);
        set(gca,'xticklabel', newlabels);

        xlabel('[a_0]');
        ylabel('\phi_1(x)');
        title('gauss1 fit (x)');
        text(0.75, 0.5, txtX1, 'FontSize', 12, 'Units', 'normalized')
        xlim(xlimit);
        ylim(ylimit);

        subplot(142)
        plot(fy, xb, yb) % data + fit

        % modify labels for tick marks
        xticks = get(gca,'xtick');
        newlabels = arrayfun(@(x) sprintf('%.1f', x), xticks, 'uni', 0);
        set(gca,'xticklabel', newlabels);

        xlabel('[a_0]');
        ylabel('\phi_1(y)');
        title('gauss1 fit (y)');
        text(0.75, 0.5, txtY1, 'FontSize', 12, 'Units', 'normalized')
        xlim(xlimit);
        ylim(ylimit);

        subplot(143)
        % X
        p31 = scatter(xa, ya, '.', 'DisplayName', 'data'); % original data
        hold on;
        p32 = plot(xdataFine, fitGausx, 'DisplayName', 'fitted curve'); % fit data
        xlabel('x array [a_0]');
        ylabel('\phi(x)');
        title('fit (x)');
        text(0.75, 0.5, txtX, 'FontSize', 12, 'Units', 'normalized')
        xlim(xlimit);
        ylim(ylimit);
        legend([p31,p32])
        
        subplot(144)
        % Y
        p41 = scatter(xb, yb, '.', 'DisplayName', 'data'); % original data
        hold on;
        p42 = plot(ydataFine, fitGausy, 'DisplayName', 'fitted curve'); % fit data
        xlabel('y array [a_0]');
        ylabel('\phi(y)');
        title('fit (y)');
        text(0.75, 0.5, txtY, 'FontSize', 12, 'Units', 'normalized')
        xlim(xlimit);
        ylim(ylimit);
        legend([p41,p42])

    end
end