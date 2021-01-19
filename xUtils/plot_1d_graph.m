function [plt] = plot_1d_graph(xlinarray,phi,gphi,tfphi,xlimit,ylimit,type,Redge,S,w,dx,datestr,delta,dims)
    plt = plot(xlinarray,phi);
    xlim([-xlimit xlimit]); ylim([0 ylimit]);
    if strcmp(type, 'x')
        xlabel('x'); ylabel('|\phi|_x^2');
    elseif strcmp(type, 'y')
        xlabel('y'); ylabel('|\phi|_y^2');
    elseif strcmp(type, 'z')
        xlabel('z'); ylabel('|\phi|_z^2');
    else
        error('Invalid type chosen. Choose from x, y, z.')
    end
    
    if nargin < 14
        dims = 2;
        warning('Number of dimensions not given. Set dims=2')
    end
    
    hold on;
    plot(xlinarray,tfphi,'LineStyle','--')
    plot(xlinarray,gphi)
    line([Redge Redge], [0 ylimit], 'Color','red','LineStyle','--'); line([-Redge -Redge], [0 ylimit], 'Color','red','LineStyle','--');
    
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
