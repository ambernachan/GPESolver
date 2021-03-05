function plot_Phis_1d(xarray, PHIS, params)
    % PHIS must be a cell array containing PhiData objects
    
    if ~iscell(PHIS)
        PHIS = {PHIS};
    end
    
    dims = length(size(PHIS{1}.phisq));
    center = size(PHIS{1}.phisq, 1) / 2;
    hold on;
    
    for i=1:length(PHIS)
        if dims == 3
            phisq = PHIS{i}.phisq;
            plot(xarray, phisq(:,center,center));
        elseif dims == 2
            phisq = PHIS{i}.phisq;
            plot(xarray, phisq(:,center));
        elseif dims == 1
            plot(xarray, PHIS{i}.phisq);
        else
            error('Number of dimensions couldn"t be determined.')
        end
    end
    
    if isfield(params, 'xlim')
        xlim([-params.xlim params.xlim]);
    end
    
    if isfield(params, 'ylim')
        ylim([0 params.ylim]);
    else
        ymax = 0;
        for i=1:length(PHIS)
            if dims == 3
                ymax = max(ymax, max(max(max(PHIS{i}.phisq))));
            elseif dims == 2
                ymax = max(ymax, max(max(PHIS{i}.phisq)));
            elseif dims == 1
                ymax = max(ymax, max(PHIS{i}.phisq));
            end
        end
        ylim([0 ymax*1.1]);
    end
    
    xlabel('simulation box axis');
    ylabel('|\phi|^2');
    
%     hold on;
%     plot(xlinarray,tfphi,'LineStyle','--')
%     plot(xlinarray,gphi)
%     line([Redge Redge], [0 ylimit], 'Color','red','LineStyle','--'); line([-Redge -Redge], [0 ylimit], 'Color','red','LineStyle','--');
%     
    % create list of chis
    chis = [];
    for i=1:length(PHIS)
        chi = PHIS{i}.S;
        chis = [chis chi];
    end
    
    dx = (max(xarray)-min(xarray)) / (size(PHIS{1}.phisq, 1) - 1);
    boxN = max(xarray);
    
    minN = Inf; maxN = 0;
    for i=1:length(PHIS)
        if dims == 3
            minN = min(minN, min(min(min(PHIS{i}.phisq))));
            maxN = max(maxN, max(max(max(PHIS{i}.phisq))));
        elseif dims == 2
            minN = min(minN, min(min(PHIS{i}.phisq)));
            maxN = max(maxN, max(max(PHIS{i}.phisq)));
        elseif dims == 1
            minN = min(minN, min(PHIS{i}.phisq));
            maxN = max(maxN, max(PHIS{i}.phisq));
        end
    end
    if maxN == minN
    
    chistr = ['S=' format_str(S) ];
    timestr = [datestr(now, 'yyyy-mm-dd') '@' datestr(now, 'HH.MM.SS')];
    boxstr = [sprintf('[-%.0g;%.0g]^%d', boxN, boxN, dims)];
    Nstr = 
    
    lgd = legend(['result S=' return_stringS(S) sprintf(',\nN=%d^%d, grid(%d,%d)^%d',...
        length(xlinarray),dims, min(xlinarray),max(xlinarray), dims)], sprintf('exp TF_{%dd} cutout',dims),...
        sprintf('exp. Gaussian\n(\\sigma\\sim%.3f,A\\sim%.2f)',w/sqrt(2),max(gphi)),...
        sprintf('R_{edge}^{%dd} = %.5f\n\t(dx = %.4f)',dims,Redge,dx));
    
    
    annotation('textbox', [0, 0.05, 0, 0], 'string', sprintf('%s',datestr))
    
    if delta == 0.5
        deltastr = '1/2';
    elseif delta >= 10
        deltastr = sprintf('%d',delta);
    elseif delta < 10 && delta >= 1
        deltastr = sprintf('%.1f',delta);
    elseif delta < 1
        deltastr = sprintf('%.3f',delta);
    else
        error('Invalid delta entered.')
    end
    
    annotation('textbox', [0, 0.95, 0, 0], 'string', ['\delta=' deltastr], 'Interpreter', 'tex', 'FitBoxToText', 'on')
end