% Plot magnetization, M = [{Mx}, {My}, {Mz}] = <Fx>, <Fy>, <Fz>
% From F = Outputs.User_defined_global(8:10)
function [] = plot_magnetizations(its, M, info, evo, method)

    close all;
    datestring = info.creationTimeString;
    if ~exist('evo','var') || isempty(evo)
        evo = 1;
    end
    
    % Determine whether this was a ground state or dynamic simulation
    % if 'ground' exists in the path name; OR neither 'ground' nor 'dynamic' does
    if ~isempty(strfind(info.name, 'ground')) || ...
            (isempty(strfind(info.name, 'ground')) && isempty(strfind(info.name, 'dynamic')))
        labelx = 'iterations';
        x = evo:evo:its*evo;
    end
    
    % if 'ground' doesn't exist in the path name but 'dynamic' does
    if isempty(strfind(info.name, 'ground'))
        if ~isempty(strfind(info.name, 'dynamic'))
            labelx = 'time (\omega_{osc}^{-1})';
            if exist('method', 'var')
                dt = real(method.Deltat);
                if imag(method.Deltat) ~= 0
%                     imag_str = ['\tau=%.2g'];
                    imag_str = sprintf('τ=%.2gΔτ', -imag(method.Deltat)/abs(method.Deltat));
                    annotation('textbox', [0.151,0.75,0.120,0.063], 'string', ...
                        sprintf('%s', imag_str), 'FitBoxToText', 'on')
                end
            else
                dt = 1;
            end
            x = evo*dt:evo*dt:its*evo*dt; % time in 1/w_ho
            [x, labelx] = makedynamictimeline(x, info);
        end
    end
    
    totalF = M{1}.^2+M{2}.^2+M{3}.^2;
    
    % marker setting step size
    numberOfMarkers = 25;
    evomarker = floor(length(x)/numberOfMarkers);
    % Creating figure
    yline(1, '--');
    hold on
    yline(-1, '--');
    plot(x, M{1}, '->', 'MarkerIndices', 1:evomarker:length(M{1}), ...
        'LineWidth', 1, 'MarkerSize', 3, 'Color', [0.6 0 0.298])
    plot(x, M{2}, '-<', 'MarkerIndices', 1:evomarker:length(M{2}), ...
        'LineWidth', 1, 'MarkerSize', 3, 'Color', [0.4 0 0.8])
    plot(x, M{3}, '-o', 'MarkerIndices', ceil(evomarker/2):evomarker:length(M{3}), ...
        'LineWidth', 1, 'MarkerSize', 3, 'Color', [0 0.4 0.8])
    plot(x, totalF, '-h', 'MarkerIndices', ceil(evomarker/2):evomarker:length(totalF), ...
        'LineWidth', 1.3, 'MarkerSize', 6, 'Color', [1 0.5 0])
    
    % Add axes labels and figure text
    xlabel(labelx); 
    ylabel('Magnetization \langleF_\alpha\rangle = \langle\Psi|F_\alpha|\Psi\rangle / |\Psi|^2');
    title('Magnetization \langleF_\alpha \rangle, and total \langleF\rangle^2');
    
    lowerlimit = max(min(min(min(min(M{1}), min(M{2})), min(M{3})), min(totalF)), -1);
    upperlimit = min(max(max(max(max(M{1}), max(M{2})), max(M{3})), max(totalF)), 1);
    if lowerlimit >= 0
        lowerlimit = -0.1;
    elseif lowerlimit < 0
        lowerlimit = lowerlimit*1.2;
    end
    if upperlimit <= 0
        upperlimit = 0.1;
    elseif upperlimit > 0
        upperlimit = upperlimit*1.2;
    end
    
    ylim([lowerlimit upperlimit]);
    
    % Ignore the constant lines in the legend
    constantLineHandles = findobj('Type', 'ConstantLine');
    for i = 1:length(constantLineHandles)
        constantLineHandles(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    
    % Add legend
    lgd = legend('\langleF_x \rangle', '\langleF_y \rangle', '\langleF_z \rangle', '\langleF\rangle^2');
    
    %add datestring to figure
    annotation('textbox', [0, 0.05, 0, 0], 'string', sprintf('%s', datestring))
    
    savename = 'Magnetization Fx,Fy,Fz over time';
    
    %add atom type text to figure
    atom_str = getAtomStr(info.params.atom);
    
%     annotation('textbox', [0.225, 0.2, 0.1, 0.1], ...
%         'string', sprintf('%s', atom_str), 'FitBoxToText', 'on')
    annotation('textbox', [0.4, 0.2, 0.1, 0.1], ...
        'string', sprintf('%s', atom_str), 'FitBoxToText', 'on', ...
        'BackgroundColor', [1, 1, 1], 'FaceAlpha', 0.8)
    
    % Save figure
    infostr = [info.fulldir(end-10:end-9) info.fulldir(end-8:end)];
    info.save_figure(1, savename, infostr)
    info.save_figure(1, savename, infostr, info.fulldir, '.png')
    hold off
 
end