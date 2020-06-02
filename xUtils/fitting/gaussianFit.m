function [A, B, U, S, uncertainties, xpeak, ypeak, f] = gaussianFit(xdata, ydata, c0, d0)
    %--------------------Gaussian Fit----------------------------------------
    %Input is xdata, ydata and must be in column  vector format
    %only two initial guesses required as the following are determine by ydata
    b0 = sqrt(max(ydata(:))-min(ydata(:)));
    a0 = min(ydata(:));
    %b0 is the guess at x-location of peak
    %c0 is guess at width of peak
    %Output gives the fitted parameters as well as xpeak location and its 
    %yvalue (i.e. xpeak & ypeak
    %-------------------------------------------------------------------------
    %Define Gauss Equation (remember the . notation
    gaussEqn ='a + b^2 .* exp(- ((x - c).^2 / (d^2) )) ';
    %Use startpoint to define guess vlaues (otherwise fit doesn't perform well)
    f = fit(xdata,ydata,gaussEqn,'Normalize','off', 'Robust' , 'bisquare', 'StartPoint',[a0, b0, c0, d0]);  %use c^2 instrad of c to enforce postive sigma/fwhm
    %f = fit(xdata,ydata,gaussEqn, 'StartPoint',[A0, B0, U0, sqrt(W0)])
    %Get x location of max value
    coeffs = coeffvalues(f);
    A = coeffs(1);
    B = coeffs(2)^2;
    U = coeffs(3);
    S = coeffs(4);
    %S = S^2;
    uncertainties = confint(f);
    db3 = abs(coeffs(2)-uncertainties(3));
    db4 = abs(coeffs(2)-uncertainties(3));
    uncertainties(3) = B - db3*2*coeffs(2); %lower B
    uncertainties(4) = B + db4*2*coeffs(2); %upper B

    % %Increase resolution of x data (by 30)
    xdataFine=(linspace(xdata(1),xdata(end),30))';  
    % Create high res Gauss function
    fitGaus =  A + B .* exp(- ((xdataFine-U).^2 ./ (S^2) ));
    %Find max value and its x location
    ypeak = max(fitGaus);
    xpeak = U;
    xpeak = mean(xpeak(:));   %Take mean incase more than one peak found
end