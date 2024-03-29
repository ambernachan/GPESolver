%% SPLitting for the Dynamic System (SPL-DS) : Explicit FFT scheme for the computation of solutions of the dynamic Gross Pitaevskii equation
%% INPUTS:
%%          Phi_0: Initial wave functions (cell array)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var1d.m)
%%          Geometry1D: Structure containing variables concerning the geometry of the problem in 1D (structure) (see Geometry1D_Var1d.m)
%%          Physics1D: Structure containing variables concerning the physics of the problem in 1D (structure) (see Physics1D_Var1d.m)
%%          Outputs: Different outputs computed during the computation of the problem (structure) (see OutputsINI_Var1D.m and FFTOutputs_Var1d.m)
%%          Print: Structure containing variables concerning the printing and drawing of informations during the computations (structure) (see Print_Var1d.m)
%%          Figure: Structure containing variables concerning the figures (structure)
%% OUTPUTS:
%%          Phi_final: Dynamic's wave functions computated with the BESP method (cell array)
%%          Outputs: Different outputs computed during the computation of the problem (structure) (see OutputsINI_Var2D.m and FFTOutputs_Var1d.m)
%% FUNCTIONS USED:
%%          FFTGeometry1D_Var1d: To compute the 1d geometry used for the FFT (line 30)
%%          FFTPhysics1D_Var1d: To compute the 1d physics used for the FFT (line 31)
%%          FFTOperators1D_Var1d: To compute the 1d operators used for the FFT (line 39)
%%          FFTOutputs_Var1d: To compute output variables (line 75,79 and 101)
%%          Print_Info1d: To print informations concerning the wave functions (line 82 and 102)
%%          Draw_solution1d: To draw the ground states (line 86)

function [Phi_final, Outputs] = SPL_DS1d(Phi_0, Method, Geometry1D, Physics1D, Outputs, Print, Figure)

%% CPUtime initilization
Method.Cputime = 0; % Initialization of the CPUtime variable
Method.Cputime_temp = cputime; % Initialization of the relative CPUtime variable

%% Geometry, operators, physics and ground states initialization for FFT
FFTGeometry1D = FFTGeometry1D_Var1d(Geometry1D); % Changing the geometry for the FFT
FFTOperators1D = FFTOperators1D_Var1d(FFTGeometry1D); % Computing the derivative FFT operators
% Changing the solution for the FFT
% FOR each component
for n = 1:Method.Ncomponents
    FFTPhi{n} = Phi_0{n}(1:FFTGeometry1D.Nx); % Removing boundaries for the FFT
end
FFTPhysics1D = FFTPhysics1D_Var1d(FFTPhi, Method, Physics1D, FFTGeometry1D, FFTOperators1D); % Changing the physics for the FFT

%% Initialization of the exponential operator for the differential operators
MatDiffx = zeros(Method.Ncomponents,FFTGeometry1D.Nx*Method.Ncomponents); % Initializing the variable that will contain the wave function with the operators applied
% FOR each component
for n = 1:Method.Ncomponents
    % FOR each component where the gradient in the x direction is non null
    for m = FFTPhysics1D.Gradientx_function_Index{n}
        MatDiffx(n,[1+(m-1)*FFTGeometry1D.Nx:m*FFTGeometry1D.Nx]) = MatDiffx(n,[1+(m-1)*FFTGeometry1D.Nx:m*FFTGeometry1D.Nx]) + FFTPhysics1D.Gradientx{n,m}.*FFTOperators1D.Gx'; % Computing the gradient operators in the x direction
    end
    MatDiffx(n,[1+(n-1)*FFTGeometry1D.Nx:n*FFTGeometry1D.Nx]) = MatDiffx(n,[1+(n-1)*FFTGeometry1D.Nx:n*FFTGeometry1D.Nx]) - FFTPhysics1D.Delta*FFTOperators1D.Dx'; % Computing the laplacian
end
FFTOperators1D.MatDiffx = MatDiffx; % Storing the operators

%% CPUtime initilization bis 
Method.Cputime_temp = cputime - Method.Cputime_temp; % Storing the CPUtime relatively to the begining of the program
Method.Cputime  = Method.Cputime + Method.Cputime_temp;  % Updating the CPUtime variable

%% Computation of dynamic of the GPE via BESP
% Stopping criterions: stopping time
while (Method.Iterations*Method.Deltat<Method.Stop_time)
    %% Updating variables
    if (Method.projection == true) && (imag(Method.Deltat) ~= 0)
        FFTPhi_tmp = FFTPhi; % Storing a temporary variable of the ground states to compute local evolution (only used for spinor BEC with imaginary time damping term!!)
    end
    Method.Cputime_temp = cputime; % Reinitialization of the relative CPUtime variable
    Method.Iterations = Method.Iterations + 1; % Iteration count

    %% Computing the dynamic system on a single time step using the splitting scheme 
    FFTPhi = Local_Splitting_solution1d(FFTPhi, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D); % Computation of the ground state using the full BESP-CNFG method on a single step of time
    
    %% Projection of the wavefunction
    % [Added by Amber de Bruijn, 2022-01-21]
    % Normalizing and ensuring conservation of longitudinal magnetization;
    % only for spinor BECs and where there is an imaginary term added to 
    % the time step, so it is not fully real.
    if Method.projection == true
        if imag(Method.Deltat) ~= 0
            if strcmp(Method.Normalization,'Multi') % normalizing the total wavefunction, no conservation of components
                Global_L2norm = 0;
                for n = 1:Method.Ncomponents
                    Global_L2norm = Global_L2norm + L2_norm1d(FFTPhi{n},FFTGeometry1D)^2; % Computing the norm of each wave function
                end

                projection = ones(3,1); projection = num2cell(projection);
                if isfield(Method, 'M') && isfield(Method, 'projection')
                    M = Method.M;    
                    if Method.projection == true
                        for n = 1:Method.Ncomponents
            %                 FFTPhi{n} = FFTPhi{n}/sqrt(Global_L2norm)*sqrt(Method.NParticles(n)); % Normalization of each wave function
                            phi{n} = sqrt(sum(abs(FFTPhi{n}).^2, 'all'));
                        end
                        pref = 1;
                        if isfield(Method, 'q') % making the scheme more numerically stable against p >> q
                            dt = imag(Method.Deltat);
                            pref = exp(-4*Method.q*dt);
                            if mod(Method.Iterations, 1000) == 0
                                sprintf('Input of prefactor in projection constant: e^(-4q*dt)=%.6g', pref)
                            end
                        end
                        projection{2} = sqrt(1 - M^2) / (sqrt( phi{2}^2 + sqrt( 4*pref*(1-M^2)*phi{1}^2*phi{3}^2 + M^2 * phi{2}^4 ) ));
                        projection{1} = sqrt( 1 + M - (projection{2}^2) * phi{2}^2 ) / (sqrt(2) * phi{1});
                        projection{3} = sqrt( 1 - M - (projection{2}^2) * phi{2}^2 ) / (sqrt(2) * phi{3});
                    end
                else % projection or M are not field of Method
                    Method.projection = false;
                end

                % projecting and calculating the global norm
                Global_L2norm = 0;
                for n = 1:Method.Ncomponents
                    FFTPhi{n} = projection{n} .* FFTPhi{n};
                    Global_L2norm = Global_L2norm + L2_norm1d(FFTPhi{n},FFTGeometry1D)^2; % Computing the norm of each wave function
                end

                % Computing the local evolution of the wavefunction, normalizing if
                % no projection constants are used.
                for n = 1:Method.Ncomponents
                    if Method.projection == false
                        FFTPhi{n} = FFTPhi{n}/sqrt(Global_L2norm)*sqrt(Method.NParticles(n)); % Normalization of each wave function
                    end
                    Method.LocalEvol(n) = max(abs(FFTPhi{n}-FFTPhi_tmp{n})); % Computing the local evolution of each wave function
                end

            elseif strcmp(Method.Normalization,'Single')

                for n = 1:Method.Ncomponents
                    FFTPhi{n} = FFTPhi{n}/L2_norm1d(FFTPhi{n},FFTGeometry1D)*sqrt(Method.NParticles(n)); % Normalization of each wave function
                    Method.LocalEvol(n) = max(abs(FFTPhi{n}-FFTPhi_tmp{n})); % Computing the local evolution of each wave function
                end
            end
        end
    end
    
    %% Updating CPUtime
    Method.Cputime_temp = cputime - Method.Cputime_temp; % Storing the CPUtime relatively to the begining of the "while" loop
    Method.Cputime = Method.Cputime + Method.Cputime_temp; % Updating the CPUtime variable
    
    %% Computation of outputs
    % IF one wants to either compute outputs or to print informations during 
    % the computation
    if (Method.Output) && (mod(Method.Iterations,Outputs.Evo_outputs) == 0)
        Outputs = FFTOutputs_Var1d(FFTPhi, Outputs, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D); % Computing the output variables
    end
    %% Printing informations and drawing ground states
    % IF the number of iterations has been reached
    if (mod(Method.Iterations,Print.Evo) == 0)
        % IF one has chosen to print informations
        if (Method.Output) && (Print.Print == 1)
            Print_Info1d(Outputs, Method) % Printing informations
        elseif (~Method.Output) && (Print.Print == 1)
            Outputstmp = FFTOutputs_Var1d(FFTPhi, Outputs, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D); % Computing the output variables
            Print_Info1d(Outputstmp, Method) % Printing informations
        end
        % IF one has chosen to draw the ground states
        if (Print.Draw == 1)
            Draw_solution1d(FFTPhi,Method,FFTGeometry1D,Figure) % Drawing the wave functions' modulus square and angle
        end;
    end
end

%% Printing the informations of the final ground states
Outputs = FFTOutputs_Var1d(FFTPhi, Outputs, Method, FFTGeometry1D, FFTPhysics1D, FFTOperators1D); % Computing the output variables
Print_Info1d(Outputs, Method) % Printing informations

%% Finalization of the ground states
% FOR each component
for n = 1:Method.Ncomponents
Phi_final{n} = zeros(1,Geometry1D.Nx); % Initializing the variable for the storing of the final ground states
Phi_final{n}(1,1:FFTGeometry1D.Nx) = FFTPhi{n}; % Storing of the ground states solutions
Phi_final{n}(FFTGeometry1D.Nx+1) = Phi_final{n}(1); % Setting the periodic boundary of the ground states solutions
end