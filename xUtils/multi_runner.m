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
    % params = Parameters(parameters);
    
    % condition for multirun
    if ~isempty(loopnames)
        for loop = 1:numel(loopnames)
            for i = 1:numel(parameters.(loopnames{loop}))
                iter = i*loop; creationTime = now;
                changedpar = parameters.(loopnames{loop})(i);
                params = parameters;
                params.(loopnames{loop}) = changedpar;
                params = Parameters(params);
                info{iter} = Info(params.scriptname, creationTime, params);
                
                % printing information about the simulation
                creationTimeString = [datestr(creationTime, 'yyyy-mm-dd') '@' datestr(creationTime, 'HH.MM.SS') ];
                sprintf('Running dynamic simulation at %s \n', creationTimeString)
                sprintf('Variable %s is set to %.5g \n', loopnames{loop}, changedpar)
                
                % run the simulation
                run_script(info{iter}, params)
            end
        end
    else
        creationTime = now;
        params = Parameters(parameters);
        info = Info(params.scriptname, creationTime, params);
        run_script(info, params)
    end
    
end