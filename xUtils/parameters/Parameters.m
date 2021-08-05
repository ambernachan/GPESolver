
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
        dx
        dt
        nComponents
        
        transitionfreq
        linewidth
        Ehfs
        detuning
        dipoleTrap0
        dipoleWaist_x
        dipoleWaist_y
        zRx
        zRy
    end
    methods
        % Constructor
        function obj = Parameters(inputparams)
            
            allprops = [{'scriptname'}, {'boxlimits'}, {'Ngridpts'},...
                {'chi'}, {'delta'}, {'gammas'}, {'dimensions'}, {'atom'},...
                {'run_dynamic'}, {'Phi_input'}, {'beta'}, {'a0'}, {'a2'},...
                {'N'}, {'trapfreq'}, {'atom_mass'}, {'spin_pair'}, {'aho'},...
                {'an'}, {'as'}, {'chin'}, {'chis'}, {'betan'}, {'betas'},...
                {'dx'}, {'dt'}, {'nComponents'} ...
                {'transitionfreq'}, {'linewidth'}, {'Ehfs'}, {'detuning'}, ...
                {'dipoleTrap0'}, {'dipoleWaist_x'}, {'dipoleWaist_y'}, ...
                {'zRx'}, {'zRy'}];
            
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
            
        end
        
        % Set the default properties of Parameters and return as a struct
        function default = setdefaultprops(obj)

            default.scriptname = 'spinor_GPE3D_dynamics';
            default.nComponents = 3;
            default.dimensions = 3;
            default.run_dynamic = true;
            default.beta = 1;
            
            boxlim = 8;
            gridpts = 2^6+1;
            default.boxlimits = [boxlim, boxlim, boxlim];
            default.Ngridpts = gridpts;
            
            default.chi = [];
            default.delta = 0.5;
            default.gammas = [1,1,1];
            
            default.atom = 'Na';
            default = obj.createAtomRelatedFields(default);
            
            default.dipoleWaist_x = getsimconst('dipole_waist_x');
            default.dipoleWaist_y = getsimconst('dipole_waist_y');
            default.zRx = getsimconst('zRx');
            default.zRy = getsimconst('zRy');
            
            default.dx = (2*boxlim) / (gridpts-1);
            default.dt = (default.dx)^2;
            
            % Creating a Thomas-Fermi Phi_input
            geom = Geometry3D_Var3d(-boxlim, boxlim, -boxlim, boxlim, -boxlim, boxlim, gridpts, gridpts, gridpts);
            potential = quadratic_potential3d(1,1,1,geom.X,geom.Y,geom.Z);
            for n = 1:default.nComponents
                Phi_in{n} = Thomas_Fermi3d(1,1,1, default.beta, potential);
            end
            default.Phi_input = Phi_in;
            
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
        end
    end
end
