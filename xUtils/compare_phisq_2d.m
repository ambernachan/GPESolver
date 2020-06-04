function compare_phisq_2d(fignum, phi1, phi2, Method, Geometry2D, varargin)

% COMPARE_PHISQ_2D 
%   compare_phisq_2d(figure number, phi1, phi2, Method, Geometry2D, varargin)
%   varargin: 'minus' / 'ratio'

difference_type = 'minus';

if nargin > 5
    difference_type = varargin{1};
end

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
    phi_ratio{n} = phi1_sq{n} ./ phi2_sq{n};
end

Function_Name = '|\phi(x,y)|^2';

%% Draw

start_fignums_phi2 = 2*Method.Ncomponents;
start_fignums_diff = 2*Method.Ncomponents*2;

phi_comparison = cell(1, Method.Ncomponents);
for n = 1:Method.Ncomponents
    if strcmp(difference_type, 'minus')
        phi_comparison{n} = phi_difference{n};
    else
        phi_comparison{n} = phi_ratio{n};
    end
end

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
    draw_function_2d(phi_comparison{n},Geometry2D,fig); % Drawing the square of the modulus of the wave function
    
    % Printing angle(phi)
    fig.label = start_fignums_diff + (fignum-1)+n+Method.Ncomponents; % Number of the figure
    fig.title = strcat('angle(\phi',32,'(x,y)) (component ', 32 ,num2str(n), '), [diff]'); % Storing title of the figure
    draw_function_2d(angle(phi1{n})-angle(phi2{n}),Geometry2D,fig); % Drawing the angle of the wave function

    
end