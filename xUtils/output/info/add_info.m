function add_info(fileID, dimensions, Method, Outputs, tElapsed, CPUtimeAtStart, CPUtimeNow, seperatorStr)

    fprintf(fileID, '%s\n', seperatorStr);
    fprintf(fileID, '%dd simulation, ', dimensions);
    if Method.Ncomponents == 1 
        fprintf(fileID, '%d component\n', Method.Ncomponents); 
    else
        fprintf(fileID, '%d components\n', Method.Ncomponents);
    end

    fprintf(fileID, 'Elapsed time:\t%s', print_time(tElapsed));
    fprintf(fileID, 'Elapsed CPU time:\t%8.2f\n', CPUtimeNow - CPUtimeAtStart);
    fprintf(fileID, 'Time:\t%s\n', datestr(now, 'dd mmm yy @ HH:MM:SS'));

    fprintf(fileID, 'max iterations:\t%d\n', Method.Max_iter);

    if isfield(Outputs,'iteration_vec')
        fprintf(fileID, '# iterations:\t%d\n', length(Outputs.iteration_vec));
    elseif Outputs.Evo_outputs == 1
        fprintf(fileID, '# iterations:\t~%d\n', Outputs.Iterations);
    else
        its = Outputs.Iterations - 0.5;
        evo = Outputs.Evo_outputs;
        fprintf(fileID, '# iterations:\t%d\x00B1%d\n', evo*its, ceil(0.5*evo));
    end

    fprintf(fileID, 'Deltat:\t%8.2e\n', Method.Deltat);
    fprintf(fileID, 'Stop time:\t%8.2f\n', Method.Stop_time);

    fprintf(fileID, '%s\n', seperatorStr);

    for i = 1:Method.Ncomponents
        fprintf(fileID, '------Outputs of component %d---------------\n', i);
        fprintf(fileID, 'Square at the origin:\t%8.14f\n', Outputs.phi_abs_0{i}(end));
        fprintf(fileID, 'x-radius mean square:\t%8.14f\n', Outputs.x_rms{i}(end));
        if dimensions > 1
            fprintf(fileID, 'y-radius mean square:\t%8.14f\n', Outputs.y_rms{i}(end));
        end
        if dimensions > 2
            fprintf(fileID, 'z-radius mean square:\t%8.14f\n', Outputs.z_rms{i}(end));
        end
        fprintf(fileID, 'Energy:\t\t\t%8.14f\n', Outputs.Energy{i}(end));
        fprintf(fileID, 'Chemical potential:\t%8.14f\n', Outputs.Chemical_potential{i}(end));
        
        if dimensions > 1
            fprintf(fileID, 'Angular momentum:\t%8.14f\n', Outputs.Angular_momentum{i}(end));
        end
        
        if strcmp(Method.Computation,'Ground')
            fprintf(fileID, 'Stopping criterion (stops at %g):\t%g\n', ...
                Method.Stop_crit{2} * Method.Deltat, Method.EvolutionCriterion);
        end
        
        % IF there are user defined local functions
        if (Outputs.User_compute_local)
        % FOR each user defined function
            for j = 1:Outputs.User_defined_number_local
                fprintf(fileID, strcat(Outputs.User_defined_names_local{j}, 32, ...
                    ': %8.14f\n'), Outputs.User_defined_local{i,j}(Outputs.Iterations));
            end 
        end
        
        % IF the number of iterations is superior to 1
        if (Outputs.Iterations>1)
            % Computing the energy decay of the wave function
            Energy_decay = Outputs.Energy{i}(Outputs.Iterations) - Outputs.Energy{i}(Outputs.Iterations-1);
        % ELSE the number of iterations is 0
        else
            Energy_decay = 0; % Setting the energy decay to zero
        end
        
        fprintf(fileID, 'Energy evolution:\t%e\n', Energy_decay);
        
    end

    fprintf(fileID, '%s\n', seperatorStr);

    if (Outputs.User_compute_global)
        % FOR each user defined function
        for j = 1:Outputs.User_defined_number_global
            fprintf(fileID, strcat(Outputs.User_defined_names_global{j}, 32, ...
                ': %8.14f\n'), Outputs.User_defined_global{j}(Outputs.Iterations));
        end
        fprintf(fileID, '%s\n', seperatorStr);
    end
end