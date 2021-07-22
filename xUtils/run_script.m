function [] = run_script(info, params)
    
    if ~params.run_dynamic
        feval(params.scriptname, info, params)
    else
        if isprop(params, 'Phi_input')
            feval(params.scriptname, info, params)
            return;
        end
        prompt = ['If you would like to change the input (ground state) function,\n',...
            'please give it as input now. Otherwise, return.\n'];
        userinput = input(prompt);
        if isempty(userinput) % leaves the Phi_input as is
            feval(params.scriptname, info, params)
        else
            params.Phi_input = userinput;
            feval(params.scriptname, info, params)
        end
    end
    
end