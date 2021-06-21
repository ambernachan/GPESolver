function params = makeparams (dimensions, boxlimits, Ngridpts, rundyn, atom, chi, delta, gammas)
    
    if nargin < 5
        params.atom = 'Na';
    end
    if nargin > 5
        params.S = chi;
        params.delta = delta;
        params.gammas = gammas;
    end
    params.dimensions = dimensions;
    params.boxlimits = boxlimits;
    params.Ngridpts = Ngridpts;
    params.dyn_simu = rundyn;
    params.atom = atom;
    
end