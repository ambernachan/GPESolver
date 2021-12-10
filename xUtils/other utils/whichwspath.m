function [wspath] = whichwspath(info)

    wspath_d = info.get_workspace_path('dynamics');
    if ~exist(wspath_d, 'file')
        wspath_g = info.get_workspace_path('groundstate');
        if exist(wspath_g, 'file')
            % case where wspath_g exists but wspath_d doesn't, return wspath_g
            wspath = wspath_g; 
        else
            % case where neither wspath_d nor wspath_g exists, so this
            % workspace was probably saved on a different computer
            laptopstr = 'D:\surfdrive\UNI\MSc Project spinor condensates\';
            computerstr = 'C:\AMBER\surfdrive\UNI\MSc Project spinor condensates\';
            
            pathfromxOut_d = [];
            p_d = split(wspath_d, [string('/'), string('\')]);
            a_d = 0;
            for i=1:length(p_d)
                if strcmp(p_d(i), 'xOutputs')
                    a_d = i;
                end
            end
            
            pathfromxOut_g = [];
            p_g = split(wspath_g, [string('/'), string('\')]);
            a_g = 0;
            for i=1:length(p_g)
                if strcmp(p_g(i), 'xOutputs')
                    a_g = i;
                end
            end
            
            if (a_d == 0) && (a_g == 0)
                % Neither wspath_d/wspath_g exist; also the string provided
                % in both cases has no 'xOutputs'. Something's wrong.
                sprintf('The atom type was not given because the workspace path was not recognized.')
                wspath = 0;
                return;
            end
               
            p_d = p_d(a_d:end);
            for i=1:length(p_d)
                pathfromxOut_d = [pathfromxOut_d '\' char(p_d(i))];
            end
            p_g = p_g(a_g:end);
            for i=1:length(p_g)
                pathfromxOut_g = [pathfromxOut_g '\' char(p_g(i))];
            end

            if strcmp(wspath_d(1), 'C')
                wspath_d = [laptopstr pathfromxOut_d];
            elseif strcmp(wspath_d(1), 'D')
                wspath_d = [computerstr pathfromxOut_d];
            end

            if strcmp(wspath_g(1), 'C')
                wspath_g = [laptopstr pathfromxOut_g];
            elseif strcmp(wspath_g(1), 'D')
                wspath_g = [computerstr pathfromxOut_g];
            end

            if exist(wspath_d, 'file')
                wspath = wspath_d;
            elseif exist(wspath_g, 'file')
                wspath = wspath_g;
            end
        end
    else
        % case where wspath_d existed
        wspath = wspath_d;
    end
end