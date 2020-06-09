function [phi, spacebox] = construct_gaussian(geometry, varargin)
    
    angle = 0; % angle (theta) given in degrees

    % Create 3d rotation matrix about the z-axis
    Rotation_matrix = @(angle) makehgtform('zrotate', angle); 
    
    if ~isfield(geometry, 'Ny') %1D
        dimensions = 1;
        [gamma_x, X0, w] = varargin{:};
        Sigma = gamma_x / (2*w^2); % Covariance matrix
        prefactor = (gamma_x / pi)^(1/4) / sqrt(w);
        %R = [1 0; 0 1]; % no rotation in 1d
    elseif ~isfield(geometry, 'Nz') %2D
        dimensions = 2;
        [gamma_x, gamma_y, X0, Y0, w] = varargin{1:5};
        if nargin > 6
            angle = varargin{6};
            if strcomp((varargin{7}),'degrees')
                angle = deg2rad(angle);
            end
        end
        R = Rotation_matrix(angle);
        R = R(1:2, 1:2); % create 2d rotation matrix
        Sigma = [gamma_x 0; 0 gamma_y]./(2*w^2);
        Sigma =  R*(Sigma)*R.';
        prefactor = (gamma_x * gamma_y)^(1/4) / (sqrt(pi)*w);
    else %3D
        dimensions = 3;
        [gamma_x, gamma_y, gamma_z, X0, Y0, Z0, w] = varargin{1:7};
        if nargin > 8
            angle = varargin{8};
            if strcomp((varargin{9}), 'degrees')
                angle = deg2rad(angle);
            end
        end
        R = Rotation_matrix(angle); % create 3d rotation matrix
        Sigma = [gamma_x 0 0; 0 gamma_y 0; 0 0 gamma_z]./(2*w^2);
        Sigma =  R*(Sigma)*R.';
        prefactor = (gamma_x * gamma_y * gamma_z)^(1/4) / ((sqrt(pi)*w)^(3/2));
    end

    %% Omega? Rotations?
    %omega = Physics.Omega;
    
    %%
    xExp = exp( -( Sigma(1)*(geometry.X-X0).^2 ) );
    
    if dimensions == 1
        phi = prefactor * xExp;
        phi = phi ./ L2_norm1d(phi, geometry); % Normalizing the Gaussian
        spacebox = {geometry.X};
    else %dim=2,3
        yExp = exp( -( Sigma(2,2)*(geometry.Y-Y0).^2 ) );
        if angle ~= 0 % add a mix term
                yExp = yExp .* exp( -( (Sigma(1,2)+Sigma(2,1)).*(geometry.X-X0).*(geometry.Y-Y0) ) );
        end
        
        if dimensions == 2
            phi = prefactor * xExp .* yExp;
            phi = phi ./ L2_norm2d(phi, geometry); % Normalizing the Gaussian
            spacebox = {geometry.X, geometry.Y};
        elseif dimensions == 3
            zExp = exp( -( Sigma(3,3)*(geometry.Z-Z0).^2 ) );
            % there are no mix terms here as the rotation is about the
            % z-axis
            phi = prefactor * xExp .* yExp .* zExp;
            phi = phi ./ L2_norm3d(phi, geometry); % Normalizing the Gaussian
            spacebox = {geometry.X, geometry.Y, geometry.Z};
        else %dim!=1,2,3
            error('Something went wrong, check function input.')
        end
    end

end
