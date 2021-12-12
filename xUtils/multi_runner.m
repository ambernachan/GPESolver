%%%%%%%%%%%%% multi_runner for running a specified code file %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% with given (or default) parameters %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% using the GPELab toolbox functions %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last updated Amber de Bruijn, 2021-12-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
DESCRIPTION OF THE MULTI_RUNNER FUNCTION
%} %{ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
===> Dependencies: format_str()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
===> Methods: 
  suffix(),                         get_workspace_path(name),
  add_simulation_info(Geometry),    save_figure(fignum, state, title,
  add_result_info(Method, Outputs),              path,extension),
  add_custom_info(str, vararg),     add_info_separator()
  finish_info(),
---> All *info*() methods are meant to add information about the
  simulations to the generated INFO.txt file. Struct arguments (Method,
  Geometry, Outputs) are assumed to be from the GPELab Toolbox.
---> get_workspace_path takes arg 'name' (e.g. 'groundstate' or 'dynamics')
  to return the full path to the workspace that is automatically saved
  upon finish_info(). 
---> save_figure requires arguments 'fignum' (figure number to be saved) 
  and 'state' (str) (e.g. 'groundstate' or 'dynamics'). Optional 
  arguments are 'title' (str) (title of the plot, left out if not given),
  'path' (str) (defaults to output path in obj), and 'extension' (str)
  (defaults to '.fig').
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for multi-parameter runs
function [info] = multi_runner(parameters)
    
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
%             for i = 1:numel(parameters.(loopnames{loop}))
%                 iter = i*loop; creationTime = now;
%                 changedpar = parameters.(loopnames{loop})(i);
%                 params = parameters;
%                 params.(loopnames{loop}) = changedpar;
%                 params = Parameters(params);
%                 info{iter} = Info(params.scriptname, creationTime, params);
                
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