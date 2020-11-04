function [] = plot_1d_graph(xlinarray,phi,gphi,tfphi,xlimit,ylimit,type,Redge,S,w,dx)
    plot(xlinarray,phi)
    xlim([-xlimit xlimit]); ylim([0 ylimit]);
    if strcmp(type, 'x')
        xlabel('x'); ylabel('|\phi|_x^2');
    else
        xlabel('y'); ylabel('|\phi|_y^2');
    end
    hold on;
    plot(xlinarray,tfphi,'LineStyle','--')
    plot(xlinarray,gphi)
    line([Redge Redge], [0 ylimit], 'Color','red','LineStyle','--'); line([-Redge -Redge], [0 ylimit], 'Color','red','LineStyle','--');
    
    lgd = legend(['result S=' return_stringS(S) sprintf(',\nN=%d^2, grid(%d,%d)^2',...
        length(xlinarray),min(xlinarray),max(xlinarray))], 'exp TF_{2d} cutout',...
        sprintf('exp. Gaussian\n(\\sigma\\sim%.3f,A\\sim%.2f)',w/sqrt(2),max(gphi)),...
        sprintf('R_{edge}^{2d} = %.5f\n\t(dx = %.4f)',Redge,dx));

end
