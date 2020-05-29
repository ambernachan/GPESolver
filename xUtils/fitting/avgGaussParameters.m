function parameters = avgGaussParameters(sCoeffs1, sCoeffs2)

    % From two GaussParameterStruct-s, generate averaged input for a new
    % GaussParameterStruct(B, X0, W, uncertainties BX0W, varargin), with
    % varargin{1} = arg#5 = parameter A; 
    % varargin{2} = arg#6 = confint A;

par = struct();
    
% Fields 1:4 of GaussParameterStruct are A,B,X0,W, which are necessary
% inputs for a new GaussParameterStruct.

fnames1 = fieldnames(sCoeffs1);
fnames2 = fieldnames(sCoeffs2);

for i = 1:4 % A B X0 W
    
    % Main field values
    fname1 = fnames1(i); fname1 = fname1{1}; % get the ith field name (1)
    fname2 = fnames2(i); fname2 = fname2{1}; % get the ith field name (2)
    ival1 = sCoeffs1.(fname1); % get the ith parameter value
    ival2 = sCoeffs2.(fname2); % get the ith parameter value

    % Find averaged value
    if (ival1 == 0 && ival2 ~= 0)
        avgVal = ival2;
        fprintf('One of the two %ss is zero.\n', num2str(fname1));
    elseif (ival1 ~= 0 && ival2 == 0)
        avgVal = ival1;
        fprintf('One of the two %ss is zero.\n', num2str(fname1));
    else
        avgVal = sqrt(ival1*ival2);
    end
    
    % save main field values in par struct
    par.(fname1) = avgVal;
    
    
    % Uncertainties
    uncfname = ['unc' fname1];
    iuncval1 = sCoeffs1.(uncfname);
    iuncval2 = sCoeffs2.(uncfname);
    % avg unc value
    if (iuncval1 == 0 && iuncval2 ~= 0)
        avgUncVal = iuncval2;
        fprintf('One of the two unc%ss is zero.\n', num2str(fname1));
    elseif (ival1 ~= 0 && iuncval2 == 0)
        avgUncVal = iuncval1;
        fprintf('One of the two unc%ss is zero.\n', num2str(fname1));
    else
        avgUncVal = sqrt(iuncval1*iuncval2);
    end
    
    upper = avgVal + abs(avgUncVal);
    lower = avgVal - abs(avgUncVal);
    
    par.(uncfname) = [upper lower];
    
end
    
%confintA = par.uncA;
%confintB = par.uncB;
%confintBX0W = [par.uncB par.uncX0 par.uncW]

parameters = GaussParameterStruct(par.B, par.X0, par.W, [par.uncB par.uncX0 par.uncW], par.A, par.uncA);
    
end