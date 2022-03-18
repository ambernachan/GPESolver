function [outputs] = fitdomainwalls(Phix, geomX, info)
    
    [xx, labl] = makexaxisinmeters(geomX, info);
    
    for i = 1:numel(Phix)
        Phisq{i} = abs(Phix{i}).^2;
    end

%     % B = 0.1 G
%     bounds{1}.neg = [892, 915];
%     bounds{2}.neg = [869, 924];
%     bounds{1}.pos = [1135, 1155];
%     bounds{2}.pos = [1125, 1178];

    % B = 1 G
    bounds{1}.neg = [885, 949];
    bounds{2}.neg = [832, 919];
    bounds{1}.pos = [1092, 1162];
    bounds{2}.pos = [1132, 1213];

    % B = 10 G
%     bounds{1}.neg = [882, 968];
%     bounds{2}.neg = [854, 918];
%     bounds{1}.pos = [1082, 1168];
%     bounds{2}.pos = [1118, 1196];
    
    X{1}.neg = geomX(bounds{1}.neg(1):bounds{1}.neg(2));
    X{2}.neg = geomX(bounds{2}.neg(1):bounds{2}.neg(2));
    X{1}.pos = geomX(bounds{1}.pos(1):bounds{1}.pos(2));
    X{2}.pos = geomX(bounds{2}.pos(1):bounds{2}.pos(2));

    x{1}.neg = xx(bounds{1}.neg(1):bounds{1}.neg(2));
    x{2}.neg = xx(bounds{2}.neg(1):bounds{2}.neg(2));
    x{1}.pos = xx(bounds{1}.pos(1):bounds{1}.pos(2));
    x{2}.pos = xx(bounds{2}.pos(1):bounds{2}.pos(2));

    PHI{1}.neg = Phisq{1}(bounds{1}.neg(1):bounds{1}.neg(2));
    PHI{2}.neg = Phisq{2}(bounds{2}.neg(1):bounds{2}.neg(2));
    PHI{1}.pos = Phisq{1}(bounds{1}.pos(1):bounds{1}.pos(2));
    PHI{2}.pos = Phisq{2}(bounds{2}.pos(1):bounds{2}.pos(2));

    % Fit equation with A,1 amplitude and shifting constants;
    % x0 is the domain (center) location and E the domain wall width
    eqn.pos = ' A .* ( tanh((x-x0)./E) +1 );';
    eqn.neg = '-A .* ( tanh((x-x0)./E) -1 );';

    % Plotting equation when we find the fit coefficient values
    formula.neg = @(E,A,x0,x) -A .* (tanh((x-x0)./E)-1);
    formula.pos = @(E,A,x0,x) A .* (tanh((x-x0)./E)+1);

    % Make fit type, setting the model using the equation defined above
    FIT.neg = fittype(eqn.neg, 'dependent', {'y'}, 'independent', {'x'}, 'coefficients', {'E', 'A', 'x0'});
    FIT.pos = fittype(eqn.pos, 'dependent', {'y'}, 'independent', {'x'}, 'coefficients', {'E', 'A', 'x0'});
    
    fitresult{1}.neg = fit(X{1}.neg, PHI{1}.neg', FIT.neg, 'StartPoint', [1, 0.001, -8.5]);
    fitresult{2}.neg = fit(X{2}.neg, PHI{2}.neg', FIT.pos, 'StartPoint', [0.60, 0.0012, -8.75]);
    fitresult{1}.pos = fit(X{1}.pos, PHI{1}.pos', FIT.pos, 'StartPoint', [1, 0.0015, 8.6]);
    fitresult{2}.pos = fit(X{2}.pos, PHI{2}.pos', FIT.neg, 'StartPoint', [0.72, 0.0013, 8.65]);
    
    plotfit{1}.neg = formula.neg(fitresult{1}.neg.E, fitresult{1}.neg.A, fitresult{1}.neg.x0, X{1}.neg);
    plotfit{2}.neg = formula.pos(fitresult{2}.neg.E, fitresult{2}.neg.A, fitresult{2}.neg.x0, X{2}.neg);
    plotfit{1}.pos = formula.pos(fitresult{1}.pos.E, fitresult{1}.pos.A, fitresult{1}.pos.x0, X{1}.pos);
    plotfit{2}.pos = formula.neg(fitresult{2}.pos.E, fitresult{2}.pos.A, fitresult{2}.pos.x0, X{2}.pos);
    
    for i = 1:2%numel(Phix)
%         Phisq{i} = abs(Phix).^2;
        
        rmse{i}.neg = sqrt(mean((PHI{i}.neg - plotfit{i}.neg').^2));
%         rmse{2}.neg = sqrt(mean((PHI{2}.neg - plotfit{2}.neg').^2));
        rmse{i}.pos = sqrt(mean((PHI{i}.pos - plotfit{i}.pos').^2));
%         rmse{2}.pos = sqrt(mean((PHI{2}.pos - plotfit{2}.pos').^2));

        meann{i}.neg = sqrt(mean((PHI{i}.neg).^2));
%         meann{2}.neg = sqrt(mean((PHI{2}.neg).^2));
        meann{i}.pos = sqrt(mean((PHI{i}.pos).^2));
%         meann{2}.pos = sqrt(mean((PHI{2}.pos).^2));

        nrmse{i}.neg = rmse{i}.neg / meann{i}.neg;
%         nrmse{2}.neg = rmse{2}.neg / meann{2}.neg;
        nrmse{i}.pos = rmse{i}.pos / meann{i}.pos;
%         nrmse{2}.pos = rmse{2}.pos / meann{2}.pos;
    end
    
    hold off;
    plot(xx, Phisq{1}, 'LineWidth', 1.2, 'Color', [0 0.6 0.1])
    hold on;
    plot(xx, Phisq{2}, 'LineWidth', 1.2, 'Color', [0 0.2 0.2])
    col = [0.8 0 0.4];
    plot(x{1}.neg, plotfit{1}.neg, 'LineStyle', 'None', 'Marker', '.', 'Color', col)
    plot(x{2}.pos, plotfit{2}.pos, 'LineStyle', 'None', 'Marker', '.', 'Color', col)
    plot(x{1}.pos, plotfit{1}.pos, 'LineStyle', 'None', 'Marker', '.', 'Color', col)
    plot(x{2}.neg, plotfit{2}.neg, 'LineStyle', 'None', 'Marker', '.', 'Color', col)
    
%     legend('\phi_+(x>0)', '\phi_0(x>0)', '\phi_+(x<0)', '\phi_0(x<0)', '|\phi_+|^2', '|\phi_0|^2')
    legend('|\phi_+|^2', '|\phi_0|^2')
    
    outputs.fitresult = fitresult;
    outputs.nrmse = nrmse;
    outputs.x = x;
    outputs.xx = xx;
    outputs.phisq = Phisq;
    outputs.plotfit = plotfit;
end