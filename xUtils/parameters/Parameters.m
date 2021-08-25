
classdef Parameters  < handle
    properties
        scriptname
        boxlimits
        Ngridpts
        chi
        delta
        gammas
        dimensions
        atom
        run_dynamic
        Phi_input
        beta
        a0
        a2
        N
        M
        trapfreq
        atom_mass
        spin_pair
        aho
        an
        as
        chin
        chis
        betan
        betas
        Bz
        Bmin % is nonzero if there is a magnetic field gradient
        p
        q
        dx
        dt
        nComponents
        projection
        
        transitionfreq
        linewidth
        Ehfs
        detuning
        dipoleTrap0
        dipoleWaist_x
        dipoleWaist_y
        zRx
        zRy
        xOmega
    end
    methods
        % Constructor
        function obj = Parameters(inputparams)
            
            allprops = [{'scriptname'}, {'boxlimits'}, {'Ngridpts'},...
                {'chi'}, {'delta'}, {'gammas'}, {'dimensions'}, {'atom'},...
                {'run_dynamic'}, {'Phi_input'}, {'beta'}, {'a0'}, {'a2'},...
                {'N'}, {'M'}, {'trapfreq'}, {'atom_mass'}, {'spin_pair'},...
                {'aho'}, {'an'}, {'as'}, {'chin'}, {'chis'}, {'betan'},...
                {'betas'}, {'Bz'}, {'Bmin'}, {'p'}, {'q'}, ...
                {'dx'}, {'dt'}, {'nComponents'}, {'projection'}, ...
                {'transitionfreq'}, {'linewidth'}, {'Ehfs'}, {'detuning'}, ...
                {'dipoleTrap0'}, {'dipoleWaist_x'}, {'dipoleWaist_y'}, ...
                {'zRx'}, {'zRy'}, {'xOmega'}];
            
            if isfield(inputparams, 'atom')
                params = obj.createAtomRelatedFields(inputparams);
            end
            
            % make a struct with all default properties
            defaultprops = obj.setdefaultprops();
            
            % Set all inputparams to obj
            inputfn = fieldnames(params);
            for k = 1:numel(inputfn)
                obj.(inputfn{k}) = params.(inputfn{k});
            end
            
            % Set all params that are still empty to default parameters
            for n = 1:numel(allprops)
                if isempty(obj.(allprops{n}))
                    obj.(allprops{n}) = defaultprops.(allprops{n});
                end
            end
            
            % Set Zeeman parameters p,q
            obj.setZeemanpars();
            obj.q = -obj.q; % make q negative (temporary)
            
            derivedQuantities(obj, inputparams);
            
        end
        
        % Set the default properties of Parameters and return as a struct
        function default = setdefaultprops(obj)

            default.scriptname = 'spinor_GPE3D_dynamics';
            default.nComponents = 3;
            default.dimensions = 3;
            default.run_dynamic = true;
            default.beta = 1;
            default.Bz = 0; % zero magnetic field
            default.Bmin = 0; % by default no magnetic field gradient
            default.p = 0; default.q = 0;
            
            boxlim = 8;
            gridpts = 2^6+1;
            default.boxlimits = [boxlim, boxlim, boxlim];
            default.Ngridpts = gridpts;
            default.M = 0;
            default.projection = true;
            
            default.chi = default.beta / (4*pi);
            default.delta = 0.5;
            default.gammas = [1,1,1];
            
            default.atom = 'Na';
            default = obj.createAtomRelatedFields(default);
            
            L = boxlim; N = gridpts;
            geom = Geometry3D_Var3d(-L,L, -L,L, -L,L, N,N,N);
            m.Ncomponents = default.nComponents;
            potential = quadratic_potential3d(1,1,1,geom.X,geom.Y,geom.Z);
            for j = 1:m.Ncomponents
                phi{j} = Thomas_Fermi3d(1,1,1, default.betan, potential);
            end
            % Normalizing
            default.Phi_input = normalize_global(m, geom, phi);
            
            default.dipoleWaist_x = getsimconst('dipole_waist_x');
            default.dipoleWaist_y = getsimconst('dipole_waist_y');
            default.zRx = getsimconst('zRx');
            default.zRy = getsimconst('zRy');
            
            default.dx = (2*boxlim) / (gridpts-1);
            default.dt = (default.dx)^2;
        end
        
        % If in the input params, certain variables are given [specifically:
        % box limits, gridpts, beta], the default struct will override
        % the quantities that are derived from these variables. This
        % function resets them according to the input.
        function derivedQuantities(obj, inputparams)
            inputprop = [{'boxlimits'}, {'Ngridpts'}, {'beta'}];

            for n = 1:numel(inputprop)
                % if the property is given in params
                if isfield(inputparams, inputprop{n})
                    % checking every derived property and inserting in obj
                    if strcmp(inputprop{n}, 'boxlimits') || strcmp(inputprop{n}, 'Ngridpts')
                        obj.(inputprop{n}) = inputparams.(inputprop{n});
                        if strcmp(inputprop{n}, 'boxlimits')
                            boxlim = inputparams.boxlimits;
                        elseif strcmp(inputprop{n}, 'Ngridpts')
                            gridpts = inputparams.Ngridpts;
                        end
                        obj.dx = 2*obj.boxlimits(1) / (obj.Ngridpts - 1);
                        obj.dt = obj.dx^2;
                    elseif strcmp(inputprop{n}, 'beta')
                        obj.beta = inputparams.beta;
                        obj.chi = obj.beta/(4*pi);
                    end
                end
            end
            
            if exist('boxlim', 'var')
                if numel(boxlim) > 1
                    L = boxlim(1);
                else
                    L = boxlim;
                end
            else
                L = obj.boxlimits(1);
            end
            if exist('gridpts', 'var')
                N = gridpts;
            else
                N = obj.Ngridpts;
            end
            geom = Geometry3D_Var3d(-L,L, -L,L, -L,L, N,N,N);
            m.Ncomponents = obj.nComponents;
            potential = quadratic_potential3d(1,1,1,geom.X,geom.Y,geom.Z);
            for j = 1:m.Ncomponents
                phi{j} = Thomas_Fermi3d(1,1,1, obj.betan, potential);
            end
            % Normalizing
            obj.Phi_input = normalize_global(m, geom, phi); 
        end
        
        function paramstruct = createAtomRelatedFields(obj, paramstruct)
            if isfield(paramstruct, 'atom')
                atom = paramstruct.atom;
            else
                atom = obj.atom;
            end
            
%             fn = fieldnames(paramstruct);
            props = [{'atom'}, {'atom_mass'}, {'N'}, {'trapfreq'}, {'spin_pair'},...
                {'aho'}, {'a0'}, {'a2'}, {'an'}, {'as'}, {'chin'}, {'chis'}, {'betan'},...
                {'betas'}, {'transitionfreq'}, {'linewidth'}, {'Ehfs'}, {'detuning'}, ...
                {'dipoleTrap0'}, {'xOmega'}];
            
            for n = 1:numel(props)
                if ~isfield(paramstruct, props{n})
                    paramstruct.(props{n}) = [];
                end
                % fill in chosen values
                if ~isempty(paramstruct.(props{n}))
                    ps.(props{n}) = paramstruct.(props{n});
                end
            end
            
            if ~isfield(ps, 'atom')
                ps.atom = atom;
            end
            if ~isfield(ps, 'atom_mass')
                ps.atom_mass = getsimconst(['mass_' ps.atom]); % Atom mass in kg
            end
            if ~isfield(ps, 'N')
                ps.N = getsimconst('N'); % number of particles
            end
            if ~isfield(ps, 'trapfreq')
                ps.trapfreq = getsimconst('trap_freq'); % Trap strength in Hz (symmetric for now)
            end
            if ~isfield(ps, 'spin_pair')
                ps.spin_pair = getsimconst('spin_pair'); % hyperfine spin manifold (=1)
            end
            
            ps.aho = sqrt(getphysconst('hbar') / (ps.atom_mass * ps.trapfreq));
            
            if ~isfield(ps, 'a0')
                ps.a0 = getsimconst(['a0_' atom]);
            end
            if ~isfield(ps, 'a2')
                ps.a2 = getsimconst(['a2_' atom]);
            end
            % TEMPORARY EDIT
            ps.a2 = ps.a0;
            
            ps.an = (2*ps.a2+ps.a0)/3;
            ps.as = (ps.a2-ps.a0)/3;
            ps.chin = ps.N*ps.an/ps.aho;
            ps.chis = ps.N*ps.as/ps.aho;
            ps.betan = 4*pi*ps.chin;
            ps.betas = 4*pi*ps.chis;
            
            pr = [{'transitionfreq'}, {'linewidth'}, {'Ehfs'}, ...
                {'detuning'}, {'dipoleTrap0'}];
            for k = 1:numel(pr)
                if ~isfield(ps, pr{k})
                    ps.(pr{k}) = getsimconst([pr{k} '_' atom]);
                end
            end
            
            if ~isfield(ps, 'xOmega')
                ps.xOmega = getsimconst(['Wx_' atom]);
            end
            
            for k = 1:numel(props)
                if isempty(paramstruct.(props{k}))
                    paramstruct.(props{k}) = ps.(props{k});
                end
            end
            
        end
        
        function [p, q] = getZeemanpars(obj, Bz, trapfreq, Ehfs)
            if nargin < 2
                Bz = obj.Bz;
                trapfreq = obj.trapfreq;
                Ehfs = obj.Ehfs;
            end
            if ~exist('Bz', 'var') || ~exist('trapfreq', 'var')
                error('Not enough input parameters; please provide Bz, Wmin, Ehfs')
            elseif ~exist('Ehfs', 'var')
                [p,q] = getMagneticFieldPars(Bz, trapfreq);
                return;
            end
            
            [p,q] = getMagneticFieldPars(Bz, trapfreq, Ehfs);
        end
        
        function setZeemanpars(obj)
            [p, q] = obj.getZeemanpars();
            obj.p = p;
            obj.q = q;
        end
    end
end
