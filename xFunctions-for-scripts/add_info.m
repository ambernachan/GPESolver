function add_info(fInfo, Outputs, Method, dimensions, tElapsed, CPUtimeAtStart, CPUtimeNow)

fprintf( fInfo, '-------------------------------------------\n');
fprintf( fInfo, '%dd simulation, ', dimensions);
if Method.Ncomponents == 1 
    fprintf( fInfo, '%d component\n', Method.Ncomponents); 
else
    fprintf( fInfo, '%d components\n', Method.Ncomponents);
end

fprintf( fInfo, 'Elapsed time:\t' );
fprintf( fInfo, print_time(tElapsed) );

fprintf( fInfo, 'max iterations:\t%d\n', Method.Max_iter);

if isfield(Outputs,'iteration_vec')
    fprintf( fInfo, '# iterations:\t%d\n', length(Outputs.iteration_vec));
elseif Outputs.Evo_outputs == 1
    fprintf( fInfo, '# iterations:\t~%d\n', Outputs.Iterations);
else
    its = Outputs.Iterations - 0.5;
    evo = Outputs.Evo_outputs;
    fprintf( fInfo, '# iterations:\t%d\x00B1%d\n', evo*its, ceil(0.5*evo));
end

fprintf( fInfo, '-------------------------------------------\n');

for i = 1:Method.Ncomponents
    fprintf( fInfo, '------Outputs of component %d---------------\n', i);
    fprintf( fInfo, 'Square at the origin:\t%8.14f\n', Outputs.phi_abs_0{i}(end));
    fprintf( fInfo, 'x-radius mean square:\t%8.14f\n', Outputs.x_rms{i}(end));
    if dimensions > 1
        fprintf( fInfo, 'y-radius mean square:\t%8.14f\n', Outputs.y_rms{i}(end));
    end
    if dimensions > 2
        fprintf( fInfo, 'z-radius mean square:\t%8.14f\n', Outputs.z_rms{i}(end));
    end
    fprintf( fInfo, 'Energy:\t\t\t%8.14f\n', Outputs.Energy{i}(end));
    fprintf( fInfo, 'Chemical potential:\t%8.14f\n', Outputs.Chemical_potential{i}(end));
    
    if dimensions > 1
        fprintf( fInfo, 'Angular momentum:\t%8.14f\n', Outputs.Angular_momentum{i}(end));
    end
    
    if strcmp(Method.Computation,'Ground')
        fprintf( fInfo, 'Stopping criterion (stops at %g):\t%g\n', ...
            Method.Stop_crit{2} * Method.Deltat, Method.EvolutionCriterion);
    end
    
    % IF there are user defined local functions
    if (Outputs.User_compute_local)
       % FOR each user defined function
        for j = 1:Outputs.User_defined_number_local
            fprintf( fInfo, strcat(Outputs.User_defined_names_local{j}, 32, ...
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
    
    fprintf( fInfo, 'Energy evolution:\t%e\n', Energy_decay);
    
end

fprintf( fInfo, '-------------------------------------------\n');

if (Outputs.User_compute_global)
    % FOR each user defined function
    for j = 1:Outputs.User_defined_number_global
        fprintf( fInfo, strcat(Outputs.User_defined_names_global{j}, 32, ...
            ': %8.14f\n'), Outputs.User_defined_global{j}(Outputs.Iterations));
    end
    fprintf( fInfo, '-------------------------------------------\n');
end

fprintf( fInfo, 'CPU time:\t%8.2f\n', CPUtimeNow - CPUtimeAtStart);
fprintf( fInfo, '-------------------------------------------\n');

%fprintf( fInfo, 'End: %s\n', datestr(now, 'dd mmm yy @ HH:MM:SS'));

end