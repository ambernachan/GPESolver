function [] = multi_runner(chi, boxlimits,...
    Ngridpts, delta, gammas, scriptname)
    creationTime = now;
    dimensions = 3;
    run_dynamic = false;

    for i = 1 : length(chi)
    %parfor i = 1 : length(chi)
        p = makeparams(chi(i), dimensions, boxlimits, Ngridpts, delta, gammas, run_dynamic);
        info = Info(scriptname, creationTime, p);
        feval(scriptname, info)
    end

end