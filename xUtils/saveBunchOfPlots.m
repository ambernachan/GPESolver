function saveBunchOfPlots(wspath, limits, opt_chi)

    close all;
    basestr = 'compto-gaussTF-';
    
%    % Assigns varname:value to current workspace
%     assignincurrent = @(varname,value) assignin('caller',varname,value);
    
    for i = 1:length(limits)
        limnames{i} = ['xlim' sprintf('%d', i)];
%         assignincurrent(limnames{i}, limits(i)); % saves xlims to current wspace
    end
    
    load(wspath)
    Rx = redge3d(1); Ry = redge3d(2); Rz = redge3d(3);
    
    ylimit_x = max(max(max(gphixN),max(tfphixN)),max(phixN))*1.1;
    ylimit_y = max(max(max(gphiyN),max(tfphiyN)),max(phiyN))*1.1;
    ylimit_z = max(max(max(gphizN),max(tfphizN)),max(phizN))*1.1;
    save(wspath, 'ylimit_x', 'ylimit_y', 'ylimit_z', '-append')
    
    % define unknowns
    if ~exist('S','var')
        S=[];
    end
    if ~exist('w','var')
        w=[];
    end
    
    %%
    % account for multi-valued S, in which case opt_chi is given:
    if length(S) > 1
        S = opt_chi;
    end
    
    % Create filenames for the figures
    for i = 1:length(limits)
        nx = [basestr 'X_' '-' format_str(limits(i)) ';' format_str(limits(i))];
        ny = [basestr 'Y_' '-' format_str(limits(i)) ';' format_str(limits(i))];
        nz = [basestr 'Z_' '-' format_str(limits(i)) ';' format_str(limits(i))];
        fnames{3*i-2} = nx; fnames{3*i-1} = ny; fnames{3*i} = nz; 
    end
    
    %%%%%%%%%%%%%%%%%%%%%

    for i = 1:length(limits)
        close all;
        
        xlim = limits(i);
        %%%%%%%%%%%%%%%%%%%%%
        plot_1d_compareG_TF(xx,phixN,gphixN,tfphixN,xlim,ylimit_x,'x',Rx,S,w,dx,datestring,Delta,info.params.dimensions)
%         plot_1d_graph(xx,phixN,gphixN,tfphixN,xlim,ylimit_x,'x',Rx,S,w,dx,datestring,Delta,info.params.dimensions)
        figname = fnames{3*i-2};
        info.save_figure(1, 'analyze', figname, fulldir, '.fig')
        info.save_figure(1, 'analyze', figname, fulldir, '.png')
        close all;

        %%%%%%%%%%%%%%%%%%%%%
        plot_1d_compareG_TF(yy,phiyN,gphiyN,tfphiyN,xlim,ylimit_y,'y',Ry,S,w,dy,datestring,Delta,info.params.dimensions)
        figname = fnames{3*i-1};
        info.save_figure(1, 'analyze', figname, fulldir, '.fig')
        info.save_figure(1, 'analyze', figname, fulldir, '.png')
        close all;

        %%%%%%%%%%%%%%%%%%%%%
        plot_1d_compareG_TF(zz,phizN,gphizN,tfphizN,xlim,ylimit_z,'z',Rz,S,w,dz,datestring,Delta,info.params.dimensions)
        figname = fnames{3*i};
        info.save_figure(1, 'analyze', figname, fulldir, '.fig')
        info.save_figure(1, 'analyze', figname, fulldir, '.png')
        close all;
    end
    
end
