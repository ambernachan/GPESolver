function params = makeparams (chi, dimensions, boxlimits, Ngridpts, delta, gammas, false)
    params.S = chi;
    params.dimensions = dimensions;
    params.boxlimits = boxlimits;
    params.Ngridpts = Ngridpts;
    params.delta = delta;
    params.gammas = gammas;
    params.dyn_simu = false;
end