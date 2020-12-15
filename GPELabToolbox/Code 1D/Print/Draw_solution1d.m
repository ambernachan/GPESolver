%% Draw wave functions' square of the modulus and angle
%% INPUTS:
%%          Phi: Wave functions (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          Geometry1D: Structure containing variables concerning the geometry of the problem in 1D (structure) (see Geometry1D_Var1d.m)
%%          Figure: Structure containing variables concerning the figures (structure) (see Figure_Var1d.m)
%% FUNCTIONS USED:
%%          draw_function_1d: To draw the wave function's square modulus and angle (line 16)

function Draw_solution1d(Phi, Method, Geometry1D, Figure)
% FOR each component
for n = 1:Method.Ncomponents
    %% Printing the figure of the square of the modulus of wave function
    Figure.label = n; % Number of the figure
    Figure.title = strcat('|phi(x)',32,'|^2 of component ', 32, num2str(n)); % Storing title of the figure
    draw_function_1d(abs(Phi{n}).^2,Geometry1D,Figure); % Drawing the square of the modulus of the wave function
    Figure.label = n + Method.Ncomponents; % Number of the figure
    Figure.title = strcat('angle(phi(x)',32,') of component ', 32, num2str(n)); % Storing title of the figure
    draw_function_1d(angle(Phi{n}),Geometry1D,Figure); % Drawing the square of the modulus of the wave function
end