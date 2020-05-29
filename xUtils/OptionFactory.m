function options = OptionFactory(varargin)

% OPTIONFACTORY Creates an options structure

%% Analysis of the inputs
Analyse_Var = inputParser; % Creating the parser

%% Adding optional arguments if they exist, otherwise setting default values

% Used in the Gaussian fit, options are 'manual' for script 'gaussianFit',
% 'gauss1' for Matlab built-in Gaussian fitting method, or 'both'
Analyse_Var.addOptional('fitmethod', 'both', @(x)ischar(x));
% Fitting option; a centered function has its peak at(/near) the origin
Analyse_Var.addOptional('centered', true, @(x)islogical(x));
% Plot all the output or simply generate it
Analyse_Var.addOptional('plotall', true, @(x)islogical(x));
% Print all the ouput or simply generate it
Analyse_Var.addOptional('printall', true, @(x)islogical(x));
% Fitting option; show the source/input image/data or not
Analyse_Var.addOptional('showsource', true, @(x)islogical(x));
% Fitting in 3d (more accurate) or not (less computation)
Analyse_Var.addOptional('dim3', false, @(x)islogical(x));

%% Parsing inputs and creating the options structure
Analyse_Var.parse(varargin{:}); % Analysing the inputs
options = Analyse_Var.Results;

end

