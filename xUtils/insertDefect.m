function [Phi] = insertDefect(Phi, defectComponents, geometry, frac)
    
    if ~exist('frac', 'var')
        frac = 0.0001;
%         frac = 0;
    end
    method.Ncomponents = length(Phi);
    geom.dx = geometry.dx; geom.dy = geometry.dx; geom.dz = geometry.dx;
    
    I = ['x', 'y', 'z'];
    for d = 1:3
        % width of the defect @ 1/100th of the simulation domain
        l = eval(['geometry.L' I(d)]) * 0.01;
        % discretized width
        w(d) = floor(l/geom.dx);
        % center coordinate
        xc(d) = (eval(['geometry.N' I(d)]) - 1) / 2;
        xmin(d) = xc(d) - w(d);
        xmax(d) = xc(d) + w(d);
    end
    
    for n = 1:length(defectComponents)
        Phi{defectComponents(n)}(xmin(2):xmax(2), xmin(1):xmax(1), xmin(3):xmax(3)) = frac * Phi{defectComponents(n)}(xmin(2):xmax(2), xmin(1):xmax(1), xmin(3):xmax(3));
    end
    
    Phi = normalize_global(method, geom, Phi);

end