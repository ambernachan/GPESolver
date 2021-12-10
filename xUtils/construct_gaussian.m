%% Constructs a rotated Gaussian function
%% Required arguments: Geometry (from GPELab struct)
%% Optional arguments: varargin:
%   1D: gamma_x, X0, w
%   2D: gamma_x, gamma_y, X0, Y0, w
%   3D: gamma_x, gamma_y, gamma_z, X0, Y0, Z0, w
%% Returns phi(gaussian), spacebox (1d: X; 2d: mesh X,Y; 3d: mesh X,Y,Z)
%% Dependencies: L2_norm1d, L2_norm2d, L2_norm3d, respectively, for 1,2,3D

function [phi, spacebox] = construct_gaussian(geometry, varargin)
    
    dimensions = 1;
    valid_angleunits = {'radians', 'degrees'}; % List of valid inputs for the angle units
    
    xInput = inputParser; % Creating the parser
    xInput.addRequired('geometry');
    xInput.addOptional('gamma_x', 1, @(x) isposrealscalar(x)); % Optional input 'gamma_x', default: 1
    xInput.addOptional('X0', 0, @(x) isscalar(x)); % Optional input 'X0' (Gaussian center), default: 0    
    if isfield(geometry, 'Ny') %>1D
        dimensions = dimensions + 1;
        xInput.addOptional('gamma_y', 1, @(x) isposrealscalar(x)); % Optional input 'gamma_y', default: 1
        xInput.addOptional('Y0', 0, @(x) isscalar(x)); % Optional input 'Y0' (Gaussian center), default: 0
        if isfield(geometry, 'Nz') %>2D
            dimensions = dimensions + 1;
            xInput.addOptional('gamma_z', 1, @(x) isposrealscalar(x)); % Optional input 'gamma_z', default: 1
            xInput.addOptional('Z0', 0, @(x) isscalar(x)); % Optional input 'Z0' (Gaussian center), default: 0
        end
    end
    if dimensions > 1
        xInput.addOptional('angle', 0, @(x) isscalar(x)); % Optional input 'angle' (rotation of the Gauss), default: 0 radians
        xInput.addOptional('angleunits', 'radians', @(x)ischar(validatestring(x, valid_angleunits))); % Optional input 'angleunits' (rotation of the Gauss), default: 'radians'
    end
    xInput.addOptional('w', 1, @(x) isposrealscalar(x)); % Optional input 'w' (width parameter), default: 1

    % Parsing inputs
    xInput.parse(varargin{:}); % Analysing the inputs
    Names = fieldnames(xInput.Results)';
    for i=1:numel(Names)
        Data{i} = xInput.Results.(Names{i});
    end
    cellfun(@(x,y) assignin('caller', x, y), Names, Data)
    
    if strcomp(xInput.Results.angleunits, 'degrees')
        angle = deg2rad(angle);
    end
    
    % Create 3d rotation matrix about the z-axis
    Rotation_matrix = @(angle) makehgtform('zrotate', angle); 
    
    if ~isfield(geometry, 'Ny') %1D
        Sigma = gamma_x / (2*w^2); % Covariance matrix
        prefactor = (gamma_x / pi)^(1/4) / sqrt(w);
        %R = [1 0; 0 1]; % no rotation in 1d
    elseif ~isfield(geometry, 'Nz') %2D
        R = Rotation_matrix(angle);
        R = R(1:2, 1:2); % create 2d rotation matrix
        Sigma = [gamma_x 0; 0 gamma_y]./(2*w^2);
        Sigma =  R*(Sigma)*R.';
        prefactor = (gamma_x * gamma_y)^(1/4) / (sqrt(pi)*w);
    else %3D
        dimensions = 3;
        R = Rotation_matrix(angle); % create 3d rotation matrix
        if angle == 0
            R = eye(3);
        end
        Sigma = [gamma_x 0 0; 0 gamma_y 0; 0 0 gamma_z]./(2*w^2);
        Sigma =  R.*(Sigma).*R';
        prefactor = (gamma_x * gamma_y * gamma_z)^(1/4) / ((sqrt(pi)*w)^(3/2));
    end

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
