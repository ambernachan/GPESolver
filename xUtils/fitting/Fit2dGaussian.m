function A = Fit2dGaussian(image, A0, sample)
    
    im = image;
    S = sample;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %                   Fit a 2D Gaussian Function to Data
    %
    %  Created by G. Nootz, May 2012. and modif by: Manuel A. Diaz, Jan 2016.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PURPOSE:  
    %   Fit a 2D gaussian centroid to simulated data.
    % 
    % INPUT:
    %   Data: two-dimensional array of size nxn.
    %   x0 = [Amp,x0,wx,y0,wy,theta]: Inital guess parameters.
    %   x = [Amp,x0,wx,y0,wy,theta]: simulated gaussian parameters.
    %   noise: noise in % of centroid peak value, A(1).
    %
    % OUTPUT: 
    %   Gaussian function parameters.
    %
    % NOTE:
    %   1.This routine uses Matlab's 'lsqcurvefit' function to fit.
    %   2.The initial values in x0 must be close to x in order for the fit
    %   to converge to the values of x (especially if noise is added).
    % 
    % MODIFS
    %   -Code is re-formulated into a single and simpler implementation.
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all; % windows

    %% ---Fitting Functions---
    %
    % Coeficients A convention:
    %	A = [Amplitude, x0, x-Width, y0, y-Width, Angle(in Radians)]
    %
    % X-data convention:
    %	X is of size(n,n,2) where 
    %	X(:,:,1) : x-coordinates,  
    %	X(:,:,2) : y-coordinates.
    %
    % In this numerical test we use two-dimensional fitting functions:

    % 1. 2D Gaussian function ( A requires 5 coefs ).
    g = @(A,X) A(1)*exp( -((X(:,:,1)-A(2)).^2/(2*A(3)^2) + (X(:,:,2)-A(4)).^2/(2*A(5)^2)) );

    % 2. 2D Rotated Gaussian function ( A requires 6 coefs ).
    f = @(A,X) A(1)*exp( -(...
        ( X(:,:,1)*cos(A(6))-X(:,:,2)*sin(A(6)) - A(2)*cos(A(6))+A(4)*sin(A(6)) ).^2/(2*A(3)^2) + ... 
        ( X(:,:,1)*sin(A(6))+X(:,:,2)*cos(A(6)) - A(2)*sin(A(6))-A(4)*cos(A(6)) ).^2/(2*A(5)^2) ) );


    %% ---Parameters---
    ss = size(im);
    n = ss(1);
    m = ss(2);
    InterpMethod='nearest'; % 'nearest','linear','spline','cubic'
    FitOrientation='dont';	% 'fit': fit for orientation, 'dont' fit for orientation

    %% ---Build numerical Grids---
    % Numerical Grid
    [x,y] = meshgrid(0:n,0:m); 
    X = zeros(m+1,n+1,2); 
    X(:,:,1) = x; 
    X(:,:,2) = y;
    % High Resolution Grid
    %h = 3; 
    h = 1;
    [xh,yh] = meshgrid(0:1/h:n,0:1/h:m); 
    Xh=zeros(h*m+1,h*n+1,2); 
    Xh(:,:,1)=xh; 
    Xh(:,:,2)=yh;

    %% ---Fit---
    % Define lower and upper bounds [Amp,xo,wx,yo,wy,fi]
    lb = [0,0,0,0,0,0];
    ub = [realmax('double'),n,(n)^2,m,(m)^2,pi/4];

    % Fit sample data
    switch FitOrientation
        case 'dont', [A,resnorm,res,flag,output] = lsqcurvefit(g,A0(1:5),X,S,lb(1:5),ub(1:5));
        case 'fit',  [A,resnorm,res,flag,output] = lsqcurvefit(f,A0,X,S,lb,ub);
        otherwise, error('invalid entry');
    end
    disp(output); % display summary of LSQ algorithm

    %% ---Plot Data---
    % Plot 3D Data and Fitted curve
    hf1 = figure(2); 
    %set(hf1,'Position',[1000 600 800 500]); 
    switch FitOrientation
        case 'dont', C=del2(g(A,Xh)); 
            mesh(xh,yh,g(A,Xh),C); hold on
        case 'fit',  C=del2(f(A,Xh)); 
            mesh(xh,yh,f(A,Xh),C); hold on
    end
    surface(x, y, S,'EdgeColor','none'); alpha(0.5); 
    colormap('pink'); %view(-60,20); 
    grid on; hold off

    % Plot Sample Pixels data
    hf2 = figure(3); 
    %set(hf2,'Position',[20 20 800 800]); 
    subplot(4,4,[5,6,7,9,10,11,13,14,15]); 
    imagesc(x(1,:),y(:,1),S); 
    colormap('hot');

    % Output and compare data and fitted function coefs
    text(0 + 20,m + 50,sprintf('\t Amplitude \t X-Coord \t X-Width \t Y-Coord \t Y-Width \t Angle'),'Color','black');
    %text(-n/2-5,m/2+6.2,sprintf('Set \t %1.3f \t %1.3f \t %1.3f \t %1.3f \t %1.3f \t %1.3f',A_s),'Color','blue');
    text(0 + 20,m + 20,sprintf('Fit \t %1.3f \t %1.3f \t %1.3f \t %1.3f \t %1.3f \t %1.3f',A),'Color','red');

    % Plot vertical and horizontal axis
    vx_h=x(1,:); vy_v=y(:,1);
    switch FitOrientation
        case 'fit', M=-tan(A(6));
            % generate points along _horizontal & _vertical axis
            vy_h = M*(vx_h-A(2))+A(4); 
            hPoints = interp2(x,y,S,vx_h,vy_h,InterpMethod);
            vx_v=M*(A(4)-vy_v)+A(2);
            vPoints = interp2(x,y,S,vx_v,vy_v,InterpMethod);
        case 'dont', A(6)=0; 
            % generate points along _horizontal & _vertical axis
            vy_h=A(4)*ones(size(vx_h)); 
            hPoints = interp2(x,y,S,vx_h,vy_h,InterpMethod);
            vx_v=A(2)*ones(size(vy_v)); 
            vPoints = interp2(x,y,S,vx_v,vy_v,InterpMethod);
    end

    % plot lines 
    hold on; 
    plot(A(2),A(4),'+b',vx_h,vy_h,'.r',vx_v,vy_v,'.g'); 
    hold off;

    % Plot cross sections 
    dmin = 1.1*min(S(:)); 
    xfit = xh(1,:); 
    hfit = A(1)*exp(-(xfit-A(2)).^2/(2*A(3)^2));
    dmax = 1.1*max(S(:)); 
    yfit = yh(:,1); 
    vfit = A(1)*exp(-(yfit-A(4)).^2/(2*A(5)^2));
    subplot(4,4,[1,2,3]); 
    xposh = (vx_h-A(2))/cos(A(6))+A(2); 
    plot(xposh,hPoints,'r.',xfit,hfit,'black'); 
    grid on; 
    axis([0,n,dmin,dmax]);
    subplot(4,4,[8,12,16]); 
    xposv = (vy_v-A(4))/cos(A(6))+A(4); 
    plot(vPoints,xposv,'g.',vfit,yfit,'black'); 
    grid on; 
    axis([dmin,dmax,0,m]); 
    set(gca,'YDir','reverse');

end