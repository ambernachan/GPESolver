function [info] = multi_runner(scriptname, boxlimits,...
    Ngridpts, chi, delta, gammas)
    creationTime = now;
    dimensions = 3;
    run_dynamic = false;
    atom = 'Rb';
    
    % if only the scriptname is given
    if nargin < 2
        boxlimits = [5,5,5];
        Ngridpts = 2^6+1;
    end
    
    if nargin > 3
        for i = 1 : length(chi)
            q(i) = makeparams(dimensions, boxlimits, Ngridpts, run_dynamic, atom, chi(i), delta, gammas);
        end
    else
        chi = 1;
        delta = 0.5;
        gammas = [1, 1, 1];
    end
    
    for i = 1 : length(chi)
    %parfor i = 1 : length(chi)
        if nargin > 3
            p = q(i);
        else
            S = [];
            p = makeparams(dimensions, boxlimits, Ngridpts, run_dynamic, atom, S, delta, gammas);
        end
        info = Info(scriptname, creationTime, p);
        if ~run_dynamic
            feval(scriptname, info)
        else
            phi_ground = input('Please enter the ground state function Phi.\n');
            if isempty(phi_ground)
                feval(scriptname, info)
            else
                feval(scriptname, info, phi_ground)
            end
        end
    end
end