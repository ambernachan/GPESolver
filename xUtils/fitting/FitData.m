
classdef PhiData  < handle
    properties
        phiData
        fittedphiData
        
        fitmethod

        guessParameters
        
        fitResultParams

        fittedphiData_X
        fittedphiData_Y
        fittedphiData_Z
    end
    methods
        % Constructor
        function obj = FitData(phi, fitmethod, guessParameters)

            if nargin < 2
                error('Too few arguments given, expected 2 or more.');
            end
            if ~isa(phi,'PhiData')
                error('Input is not a PhiData structure.');
            end

            obj.phiData = phi;
            obj.fitmethod = fitmethod;
            if nargin > 2
                obj.guessParameters = guessParameters;
            end

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
        
        function [] = func(obj)

            [fitresult, zfit, fiterr, zerr, resnorm, rr] = fmgaussfit(phiData.X,phiData.Y,phiData.phisq)
        end

    end
    methods (Access = private)

        function [] = func(obj)
            a = 1;
        end

    end

end
