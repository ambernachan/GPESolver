function params = makeparams (dimensions, boxlimits, Ngridpts, false, chi, delta, gammas)
    if nargin > 4
        params.S = chi;
        params.delta = delta;
        params.gammas = gammas;
    end
    params.dimensions = dimensions;
    params.boxlimits = boxlimits;
    params.Ngridpts = Ngridpts;
    params.dyn_simu = false;
end