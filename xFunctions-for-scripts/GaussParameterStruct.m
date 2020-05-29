function coeffs = GaussParameterStruct(B, X0, W, uncertaintiesBX0W, varargin)

%{
arg5; varargin{1} = parameter A
arg6; varargin{2} = confint for coefficient A
%}

% y = a + b*exp(-(x-x0)^2/w^2) fits coeffs A,B,U,W (=a, b, x0, w)

coeffs = struct();

if nargin > 4
    coeffs.A = varargin{1};
else
    coeffs.A = 0;
end

coeffs.B = B;
coeffs.X0 = X0;
coeffs.W = W;
coeffs.Sigma = W/sqrt(2); %std deviation
coeffs.Var = coeffs.Sigma^2; %variance

fnames = fieldnames(coeffs);
siz = size(fnames); siz = siz(1); % number of fields in coeffs

% Set parameters+uncertainties to zero when they have imaginary values
for i = 1:siz % loop over parameters in coeffs struct
    fname = fnames(i); fname = fname{1}; % get the ith parameter name
    ival = coeffs.(fname); % get the ith parameter value
    
    if imag(ival) ~= 0 % imaginary parameter value
        coeffs.(fname) = 0; % set parameter value to zero
        coeffs.(['unc' fname]) = 0; % create and set uncertainty=0
    end
    
    % for non-imaginary field values, set the uncertainties for ABX0WSV
    if ~isfield(coeffs,(['unc' fname])) 
        
        if i == 1 % A
            if nargin > 5
                intA = varargin{2};
                upper = abs(coeffs.A-intA(1));
                lower = abs(coeffs.A-intA(2));
            else
                upper = 0;
                lower = 0;
            end
            coeffs.uncA = max(upper, lower);
        elseif i>1 && i<5 % B X0 W
            int = uncertaintiesBX0W((i*2-3):(i*2-2));
            upper = abs(coeffs.(fname) - int(1));
            lower = abs(coeffs.(fname) - int(2));
            coeffs.(['unc' fname]) = max(upper,lower);
        elseif i == 5 % sigma
            coeffs.uncSigma = coeffs.uncW/sqrt(2);
        elseif i == 6 % variance
            coeffs.uncVar = 2 * coeffs.Sigma * coeffs.uncSigma;
        end
    end
    
    % Determine relative uncertainties
    if ival == 0 
        coeffs.(['relativeUnc' fname]) = 0;
    else
        coeffs.(['relativeUnc' fname]) = coeffs.(['unc' fname]) / ival;
    end
    
end    

% for very small X0 (probably centered Gaussians), set the relative
% uncertainty to zero, as otherwise the ratio uncX0/X0 will be huge
if abs(coeffs.X0) < 0.001*abs(B)
    coeffs.relativeUncX0 = 0;
else
    coeffs.relativeUncX0 = coeffs.uncX0/coeffs.X0;
end
    
coeffs = orderfields(coeffs, {...
    'A', 'B', 'X0', 'W', 'Sigma', 'Var', 'uncA',...
    'uncB', 'uncX0', 'uncW', 'uncSigma', 'uncVar',...
    'relativeUncA', 'relativeUncB', 'relativeUncX0',...
    'relativeUncW', 'relativeUncSigma', 'relativeUncVar'});

end