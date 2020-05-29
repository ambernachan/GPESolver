function draw_phisq_2d(fignum, phi, Method, Geometry2D, varargin)

% Create the function of phi that must be drawn
Default_Function = @(phi,X,Y) abs(phi).^2;

% Create a Figure
fig = Figure_Var2d();
fig.label = fignum; % unnecessary?
fig.title = '';

%% Analysis of inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('Function_Name','|\phi(x,y)|^2',@(x) (ischar(x)) || (iscell(x))); % Optional input 'Function_Name'
Analyse_Var.addOptional('Function', Default_Function); % Optional input 'Function'

%% Parsing inputs
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
% Setting inputs
Function = Analyse_Var.Results.Function; %Storing the 'Function' input
Function_Name = Analyse_Var.Results.Function_Name; %Storing the 'Function_Name' input

%%

for n = 1:Method.Ncomponents
    %% Printing the figure of the square of the modulus of wave function
    fig.label = (fignum-1)+n; % Number of the figure
    fig.title = strcat(Function_Name,' of component ', 32, num2str(n)); % Storing title of the figure
    draw_function_2d(Function(phi{n},Geometry2D.X,Geometry2D.Y),Geometry2D,fig); % Drawing the square of the modulus of the wave function
    %% Printing the figure of the angle of wave function
    fig.label = (fignum-1)+n+Method.Ncomponents; % Number of the figure
    fig.title = strcat('angle(\phi',32,'(x,y)) of component ', 32 ,num2str(n)); % Storing title of the figure
    draw_function_2d(angle(phi{n}),Geometry2D,fig); % Drawing the angle of the wave function
end
