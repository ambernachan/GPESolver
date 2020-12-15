%% Draw wave function
%% INPUTS:
%%          phi: Wave function (matrix)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          Geometry2D: Structure containing variables concerning the geometry of the problem in 2D (structure) (see Geometry2D_Var2d.m)
%%          Figure: Structure containing variables concerning the figures (structure) (see Figure_Var2d.m)

function draw_function_2d(phi,Geometry2D,Figure)
figure(Figure.label); % Setting the number of the figure
clf(Figure.label); % Clear figure
surf(Geometry2D.X,Geometry2D.Y,phi,'EdgeColor','none'); % Drawing function
shading interp; % Setting shading
colormap(Figure.map); % Setting colormap
colorbar; % Setting colorbar
view(2); % Setting view
xlabel(Figure.x); % Setting x-axis label
ylabel(Figure.y); % Setting y-axis label
title(Figure.title); % Setting title of the figure
drawnow; % Drawing