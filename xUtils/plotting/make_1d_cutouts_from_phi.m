function [F] = make_1d_cutouts_from_phi(phi)
    
    dims = ndims(phi); % 2d/3d
    
    %% find max columns
    if dims == 1
        fx = phi;
        F = {fx};
    elseif dims == 2
        
        fx = sum(phi,2);
        [~,col] = max(fx);
        fx = phi(col,:);
        
        fy = sum(phi,1);
        [~,col] = max(fy);
        fy = phi(:,col);
        
        F = {fx,fy};
        
    elseif dims == 3
        fx = sum(sum(phi,1),3);
        fy = sum(sum(phi,2),3);
        fz = sum(sum(phi,2),1);
        
        [~,colx] = max(fx);
        [~,coly] = max(fy);
        [~,colz] = max(fz);

        fx = phi(coly,:,colz);
        fy = phi(:,colx,colz)';
        fz = phi(colx,coly,:); fz = fz(:)';
        
        F = {fx,fy,fz};
    else
        error('dimensions not found');
    end


end