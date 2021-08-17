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