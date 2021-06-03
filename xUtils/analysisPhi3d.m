function [xx, phix, tfphix, gphix] = analysisPhi3d(wspath)

    close all;
    
    load(wspath) % contains Methods, Phis, simulation parameters & paths

    % make 1d slices out of the result phisq
    if exist('phi_dyn', 'var')
        phi = phi_dyn;
    else
        phi = phi_ground;
    end
    PHIS = make_1d_cutouts_from_phi(phi.phisq{1});
    phix = PHIS{1}; phiy = PHIS{2}; phiz = PHIS{3};
    geom = phi.return_geometry();
    
    % create expectation functions
    if ~exist('S','var')
        expphis = ExpPhis(1, [gx,gy,gz], [0,0,0], geom);
    else 
        expphis = ExpPhis(S, [gx,gy,gz], [0,0,0], geom);
    end
    
    
    [g1d,wg1d] = expphis.Gaussian1d();
    [g2d,wg2d] = expphis.Gaussian2d();
    [g3d,wg3d] = expphis.Gaussian3d();
    [tf3d,redge3d] = expphis.TF3d();
    [tf2d,redge2d] = expphis.TF2d();
    [tf1d,redge1d] = expphis.TF1d();
    
    % make 1d slices out of the expectation functions
    tfPHIS = make_1d_cutouts_from_phi(tf3d.phisq{1}); 
    tfphix = tfPHIS{1}; tfphiy = tfPHIS{2}; tfphiz = tfPHIS{3};
    gPHIS = make_1d_cutouts_from_phi(g3d.phisq{1}); 
    gphix = gPHIS{1}; gphiy = gPHIS{2}; gphiz = gPHIS{3};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create simulation box slice arrays & find resolution
    xx = geom.X(1,:,1);
    yy = geom.Y(:,1,1)';
    zz = geom.Z(1,1,:); zz = zz(:);
    dx = geom.dx;
    dy = geom.dy;
    dz = geom.dz;
    
    % normalize the slices
    areax = (dx * sum(phix));
    areay = (dy * sum(phiy));
    areaz = (dz * sum(phiz));
    
    phixN = phix / areax;
    phiyN = phiy / areay;
    phizN = phiz / areaz;

    gphixN = gphix / (dx*sum(gphix));
    gphiyN = gphiy / (dy*sum(gphiy));
    gphizN = gphiz / (dz*sum(gphiz));
    tfphixN = tfphix / (dx*sum(tfphix));
    tfphiyN = tfphiy / (dy*sum(tfphiy));
    tfphizN = tfphiz / (dz*sum(tfphiz));

    % save all results to the workspace
    save(wspath, '-append')
    
end