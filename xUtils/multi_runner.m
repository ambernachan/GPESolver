%%%%%%%%%%%%% multi_runner for running a specified code file %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% with given (or default) parameters %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% using the GPELab toolbox functions %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last updated Amber de Bruijn, 2021-12-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
DESCRIPTION OF THE MULTI_RUNNER FUNCTION
===> MULTI_RUNNER(parameters) takes a struct list of parameters with values
        and uses them as input for creating a Parameters object and an Info
        object, which make sure that all parameter info and data is stored
        in an identical way, making data tracking and analysis easier. The
        Parameters object is also initialized with default values for each
        parameter that is not explicitly given, so that there is less
        chance of forgetting parameters (may lead to errors, so first
        checking your list of parameters by trying parameters =
        Parameters(parameters) is advisable).
---> This function checks for tuples in the list of parameters, so that a
        multiple-variable run is possible, with exception variables
        "gammas", "boxlimits", "Phi_input", "scriptname" and "atom", which
        are expected to be tuples for single-runs. Other parameters, e.g.
        number of grid points "N" or others, will trigger multiple runs to
        start, each with different values of the input parameters.
%} %{ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
===> Dependencies: 
        B = allcomb(A1,A2,...,AN) returns all combinations of the elements 
            in the arrays A1, A2, ..., and AN.
        num = elementIndex(nthElement, elementNames, multiplicity,
            currentRunNumber) returns the index of the value of a to-be-
            varied element in a list of elements. With a list of element 
            names and the multiplicity of those elements (i.e. ["A1", "A2",
            "A3"] and [4, 5, 6], respectively) along with a 
            currentRunNumber (i.e. up to 120) and an nthElement (1, 2, or 
            3): num can be 1-4 for nthElement=1; 1-5 for nthElement=2; and 
            1-6 for nthElement=3.
        parameters = Parameters(parameterstruct) returns the initialized
            (with default values for various parameters necessary for 
            running the simulation script)
        info = Info(scriptname (str), creationTime (int), parameters 
            (Parameters class)) returns an Info object that contains and
            automatically saves simulation information in an info.txt file
            and Matlab workspaces.
        run_script(info (Info class), parameters (Parameters class)) runs
            the file indicated by Info info (where the script name is 
            stored, among other parameters).
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [info] = multi_runner(parameters)
    
    % Initialize variables to be able to detect tuples (when we need to do
    % multiple runs)
    fn = fieldnames(parameters);
    loopnames = [];
    nLoop = 0;
    
    % fieldnames that are allowed to be tuples
    tuples = struct('gammas', [], 'boxlimits', [], 'Phi_input', [], ...
        'scriptname', [], 'atom', []);
    for k = 1:numel(fn)
        if ~isfield(tuples, fn{k}) && length(parameters.(fn{k})) > 1
                % create a new loop
                nLoop = nLoop + 1;
                loopnames = [loopnames, {fn{k}}];
        else
            continue;
        end
    end
    
    % condition for multirun
    if ~isempty(loopnames)
        
        % Create a list of all combinations of parameters
        runlist = []; multiplicity = [];
        for k = 1:numel(loopnames)
            runlist = [runlist {parameters.(loopnames{k})}];
            multiplicity = [multiplicity {numel(parameters.(loopnames{k}))}];
        end
        multiplicity = [multiplicity{:}];
        runlist = allcomb(runlist{:}); % all combinations
        
        % loop over the combination list
        for loop = 1:length(runlist)
            sprintf('Run #%d/%d\n', loop, length(runlist))
            params = parameters;
            creationTime = now;
            for k = 1:numel(loopnames)
                elementIdx = elementIndex(k, loopnames, multiplicity, loop);
                changedpar = parameters.(loopnames{k})(elementIdx);
                params.(loopnames{k}) = changedpar;
                params = Parameters(params);
                info{loop} = Info(params.scriptname, creationTime, params);
                sprintf('Variable %s is set to %.5g \n', loopnames{k}, changedpar)
            end
                
            % printing information about the simulation
            creationTimeString = [datestr(creationTime, 'yyyy-mm-dd') '@' datestr(creationTime, 'HH.MM.SS') ];
            sprintf('Running simulation at %s \n', creationTimeString)
            
            % run the simulation
            run_script(info{loop}, params)
        end
    else
        creationTime = now;
        params = Parameters(parameters);
        info = Info(params.scriptname, creationTime, params);
        run_script(info, params)
    end
    
end