%% Draw wave functions' square of the modulus and angle
%% INPUTS:
%%          Phi: Wave functions (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          Geometry2D: Structure containing variables concerning the geometry of the problem in 2D (structure) (see Geometry2D_Var2d.m)
%%          Figure: Structure containing variables concerning the figures (structure) (see Figure_Var2d.m)
%% FUNCTIONS USED:
%%          draw_function_2d: To draw the wave function's square modulus and angle (line 19 and 23)

function Draw_solution2d(Phi, Method, Geometry2D, Figure)
% FOR each component
for n = 1:Method.Ncomponents
    %% Printing the figure of the square of the modulus of wave function
    Figure.label = n; % Number of the figure
    Figure.title = strcat('|phi(x,y)|^2 of component ', 32, num2str(n)); % Storing title of the figure
    draw_function_2d(abs(Phi{n}).^2,Geometry2D,Figure); % Drawing the square of the modulus of the wave function
    %% Printing the figure of the angle of wave function
    Figure.label = n+Method.Ncomponents; % Number of the figure
    Figure.title = strcat('angle(\phi',32,'(x,y)) of component ', 32 ,num2str(n)); % Storing title of the figure
    draw_function_2d(angle(Phi{n}),Geometry2D,Figure); % Drawing the angle of the wave function
end