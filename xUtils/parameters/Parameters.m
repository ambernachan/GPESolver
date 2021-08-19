
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
                {'betas'}, {'Bz'}, {'p'}, {'q'}, ...
                {'dx'}, {'dt'}, {'nComponents'}, {'projection'}, ...
                {'transitionfreq'}, {'linewidth'}, {'Ehfs'}, {'detuning'}, ...
                {'dipoleTrap0'}, {'dipoleWaist_x'}, {'dipoleWaist_y'}, ...
                {'zRx'}, {'zRy'}, {'xOmega'}];
            
            if isfield(inputparams, 'atom')
                inputparams = obj.createAtomRelatedFields(inputparams);
            end
            
            % make a struct with all default properties
            defaultprops = obj.setdefaultprops();
            
            inputfn = fieldnames(inputparams);
            for k = 1:numel(inputfn)
                obj.(inputfn{k}) = inputparams.(inputfn{k});
            end
            
            for n = 1:numel(allprops)
                if isempty(obj.(allprops{n}))
                    obj.(allprops{n}) = defaultprops.(allprops{n});
                end
            end
            
            obj.setZeemanpars();
            obj.q = -obj.q; % make q negative (temporary)
            
        end
        
        % Set the default properties of Parameters and return as a struct
        function default = setdefaultprops(obj)

            default.scriptname = 'spinor_GPE3D_dynamics';
            default.nComponents = 3;
            default.dimensions = 3;
            default.run_dynamic = true;
            default.beta = 1;
            default.Bz = 0; % zero magnetic field
            default.p = 0; default.q = 0;
            
            boxlim = 8;
            gridpts = 2^6+1;
            default.boxlimits = [boxlim, boxlim, boxlim];
            default.Ngridpts = gridpts;
            default.M = 0;
            default.projection = true;
            
            default.chi = [];
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
        
        function paramstruct = createAtomRelatedFields(obj, paramstruct)
            if isfield(paramstruct, 'atom')
                atom = paramstruct.atom;
            else
                atom = obj.atom;
            end
            paramstruct.atom_mass = getsimconst(['mass_' paramstruct.atom]); % Atom mass in kg
            paramstruct.N = getsimconst('N'); % number of particles
            paramstruct.trapfreq = getsimconst('trap_freq'); % Trap strength in Hz (symmetric for now)
            paramstruct.spin_pair = getsimconst('spin_pair'); % hyperfine spin manifold (=1)
            paramstruct.aho = sqrt(getphysconst('hbar') / (paramstruct.atom_mass * paramstruct.trapfreq));
            
            paramstruct.a0 = getsimconst(['a0_' atom]);
            paramstruct.a2 = getsimconst(['a2_' atom]);
            paramstruct.an = (2*paramstruct.a2+paramstruct.a0)/3;
            paramstruct.as = (paramstruct.a2-paramstruct.a0)/3;
            
            paramstruct.chin = paramstruct.N*paramstruct.an/paramstruct.aho;
            paramstruct.chis = paramstruct.N*paramstruct.as/paramstruct.aho;
            paramstruct.betan = 4*pi*paramstruct.chin;
            paramstruct.betas = 4*pi*paramstruct.chis;
            
            paramstruct.transitionfreq = getsimconst(['transitionfreq_' atom]);
            paramstruct.linewidth = getsimconst(['linewidth_' atom]);
            paramstruct.Ehfs = getsimconst(['Ehfs_' atom]);
            paramstruct.detuning = getsimconst(['detuning_' atom]);
            paramstruct.dipoleTrap0 = getsimconst(['dipoleTrap0_' atom]);
            paramstruct.xOmega = getsimconst(['Wx_' atom]);
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
