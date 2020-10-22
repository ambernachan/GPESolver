
classdef PhiData  < handle
    properties
        phiData
        fittedphiData
        
        fitmethod
        fitfunc

        guessParameters
        
        fitResultParams

        fittedphiData_X
        fittedphiData_Y
        fittedphiData_Z
    end
    methods
        % Constructor
        function obj = FitData(phi, fitmethod, guessParameters)

            % if only 'phi' is provided:
            if nargin < 2
                prompt = 'Give fitmethod: "gauss", "thomasfermi", or "function".';
                obj.fitmethod = input(prompt, 's');
                if strcmp(obj.fitmethod, 'function')
                    prompt = 'Provide the fit function.';
                    obj.fitfunc = input(prompt);
                end
                % if a function is provided rather than the input 'function':
                if ~strcmp(obj.fitmethod, 'gauss') && ~strcmp(obj.fitmethod, 'thomasfermi') && ~strcmp(obj.fitmethod, 'function')
                    obj.fitfunc = obj.fitmethod;
                    obj.fitmethod = 'function';
                end
            end

            if ~isa(phi,'PhiData')
                error('Input is not a PhiData structure.');
            end
            obj.phiData = phi;

            % if 'fitmethod' is provided
            if nargin > 1
                % If the 'fitmethod' is not given as options 'gauss' or 'thomasfermi', then it
                % must be 'function', and a function may be given in its place, so the input is 
                % saved as the obj.fitfunction and the obj.fitmethod is set to be 'function'.
                if ~strcmp(fitmethod, 'gauss') && ~strcmp(fitmethod, 'thomasfermi')
                    obj.fitfunc = fitmethod;
                    obj.fitmethod = 'function';
                end
            end
            
            % if 'guessParameters' is provided
            if nargin > 2
                obj.guessParameters = guessParameters;
            end

            % ------------------------------------------------------------------ fix below ---------

            if strcmp(fitmethod, '1dgauss')
                if nargin < 3
                    error('Too few arguments given; %s fit requires guessParameters', obj.fitmethod);
                end
            elseif strcmp(fitmethod, 'fullfit')
                if nargin > 2
                    error('guessParameters given but not used.');
                end
                if obj.phiData.dimensions == 1
                    obj.fitmethod = 'gauss1';
                end
                if obj.phiData.dimensions == 3
                    error('fitmethod fullfit not possible for 3D simulations.')
            elseif (~isstring(fitmethod) || ~ischar(fitmethod))
                error('fitmethod is not a string or char array');
            elseif ~strcmp(fitmethod, 'gauss1') % for custom fits with given function
                if nargin < 3
                    error('Too few arguments given; custom fit requires guessParameters');
                end
            end

            % Options for obj.fitmethod:
            %   'gauss1',   built-in Matlab function for fitting a 1D array
            %                  This will yield one result for 1D sim; two for 2D sim
            %                   and three for 3D sim (= xarray fit (& yarray fit (& zarray fit))),
            %                   along with an averaged result. 
            %                   Performs a 1d fit with fit function 'a*exp(-((x-b)/c)^2)'
            %                   * OPTIONAL guessParameters
            %   '1dgauss'   In a similar manner as 'gauss1', so yielding more results for higher-d
            %                   simulations, uses script 'gaussianFit', which:
            %                   Performs a 1d fit with fit function 'a + b^2 .* exp(- (x-c).^2/(d^2) )'. 
            %                   * DOES require guessParameters.
            %   'fullfit',  this fit method is accessible for phi input that is created in a 2D simulation, 
            %                   defaults to 'gauss1' for a 1D simulation and there is no solution yet for
            %                   3D simulation data. 
            %                   * DOESN'T require guessParameters
            %   'f(X,Y,{abc})', All other input will be interpreted as a fitfunction, i.e. 'a*x.^2 + b*x + c'.
            %                   * DOES require (i.e.) guessParameters = [a_0, b_0, c_0]
           
        end

        %% Fitting
        % Note that |phi|^2 is fitted!

        function fullfit(obj)
            
            [fitresult, zfit, fiterr, zerr, resnorm, rr] = fmgaussfit(phiData.X,phiData.Y,phiData.phisq)
        end

        function gauss1(obj)
            if obj.phiData.dimensions == 1
                X = cell(1, obj.phiData.ncomponents);
                xaxis = cell(1, obj.phiData.ncomponents);
                fitresult = cell(1, obj.phiData.ncomponents)
                for n = 1:obj.ncomponents
                    X{n} = obj.phiData.phisq{n};
                    xaxis{n} = obj.phiData.X;
                end
            end
            if obj.phiData.dimensions == 2
                X = cell(1, obj.phiData.ncomponents);
                Y = cell(1, obj.phiData.ncomponents);
                xaxis = cell(1, obj.phiData.ncomponents);
                yaxis = cell(1, obj.phiData.ncomponents);
                fitresult = cell(1, obj.phiData.ncomponents)
                for n = 1:obj.ncomponentsv
                    (X{n}, Y{n}) = obj.phiData.phisq{n};
                    xaxis{n} = obj.phiData.X;
                    yaxis{n} = obj.phiData.Y;
                end
            end
            if obj.phiData.dimensions == 3
                X = cell(1, obj.phiData.ncomponents);
                Y = cell(1, obj.phiData.ncomponents);
                Z = cell(1, obj.phiData.ncomponents);
                xaxis = cell(1, obj.phiData.ncomponents);
                yaxis = cell(1, obj.phiData.ncomponents);
                zaxis = cell(1, obj.phiData.ncomponents);
                fitresult = cell(1, obj.phiData.ncomponents)
                for n = 1:obj.ncomponents
                    (X{n}, Y{n}, Z{n}) = obj.phiData.phisq{n};
                    xaxis{n} = obj.phiData.X;
                    yaxis{n} = obj.phiData.Y;
                    zaxis{n} = obj.phiData.Z;
                end
            end

            sprintf('Starting fit (test message)...');

            

            f = fit(x.',y.','gauss1')
        end

    end
    methods (Access = private)

        function [] = func(obj)
            a = 1;
        end

    end

end
