function compare_phisq_2d(fignum, phi1, phi2, Method, Geometry2D, varargin)

% COMPARE_PHISQ_2D 
%   compare_phisq_2d(figure number, phi1, phi2, Method, Geometry2D, varargin)
% 
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

% Create a Figure
fig = Figure_Var2d();
fig.label = fignum; % unnecessary?
fig.title = '';

% Create phi squared (density) functions
Default_Function = @(phi,X,Y) abs(phi).^2;

phi1_sq = cell(1,Method.Ncomponents);
phi2_sq = cell(1,Method.Ncomponents);
phi_difference = cell(1,Method.Ncomponents);
phi_reverse_difference = cell(1,Method.Ncomponents);

for n = 1:Method.Ncomponents
    phi1_sq{n} = Default_Function(phi1{n}, Geometry2D.X, Geometry2D.Y);
    phi2_sq{n} = Default_Function(phi2{n}, Geometry2D.X, Geometry2D.Y);
    phi_difference{n} = phi1_sq{n} - phi2_sq{n};
    phi_reverse_difference{n} = phi2_sq{n} - phi1_sq{n};
end

Function_Name = '|\phi(x,y)|^2';

%% Draw

start_fignums_phi2 = 2*Method.Ncomponents;
start_fignums_diff = 2*Method.Ncomponents*2;

for n = 1:Method.Ncomponents
    
    %% phi_1
    % Printing phi^2
    fig.label = (fignum-1)+n; % Number of the figure
    fig.title = strcat(Function_Name,' (component ', 32, num2str(n), '), [a]'); % Storing title of the figure
    draw_function_2d(phi1_sq{n},Geometry2D,fig); % Drawing the square of the modulus of the wave function
    
    % Printing angle(phi)
    fig.label = (fignum-1)+n+Method.Ncomponents; % Number of the figure
    fig.title = strcat('angle(\phi',32,'(x,y)) (component ', 32 ,num2str(n), '), [a]'); % Storing title of the figure
    draw_function_2d(angle(phi1{n}),Geometry2D,fig); % Drawing the angle of the wave function

    %% phi_2
    % Printing phi^2
    fig.label = start_fignums_phi2 + (fignum-1)+n; % Number of the figure
    fig.title = strcat(Function_Name,' (component ', 32, num2str(n), '), [b]'); % Storing title of the figure
    draw_function_2d(phi2_sq{n},Geometry2D,fig); % Drawing the square of the modulus of the wave function
    
    % Printing angle(phi)
    fig.label = start_fignums_phi2 + (fignum-1)+n+Method.Ncomponents; % Number of the figure
    fig.title = strcat('angle(\phi',32,'(x,y)) (component ', 32 ,num2str(n), '), [b]'); % Storing title of the figure
    draw_function_2d(angle(phi2{n}),Geometry2D,fig); % Drawing the angle of the wave function
    
    %% phi_diff
    % Printing phi^2
    fig.label = start_fignums_diff + (fignum-1)+n; % Number of the figure
    fig.title = strcat(Function_Name,' (component ', 32, num2str(n), '), [diff]'); % Storing title of the figure
    draw_function_2d(phi_difference{n},Geometry2D,fig); % Drawing the square of the modulus of the wave function
    
    % Printing angle(phi)
    fig.label = start_fignums_diff + (fignum-1)+n+Method.Ncomponents; % Number of the figure
    fig.title = strcat('angle(\phi',32,'(x,y)) (component ', 32 ,num2str(n), '), [diff]'); % Storing title of the figure
    draw_function_2d(angle(phi1{n})-angle(phi2{n}),Geometry2D,fig); % Drawing the angle of the wave function

    
end