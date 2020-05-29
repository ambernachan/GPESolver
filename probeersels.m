clear;

%% CHOOSE DIMENSIONALITY TO TEST
dims = 1;

%% Saving info files

workingDir = pwd;
workingFile = mfilename;
fnameBase = workingFile(1:11); % = 'probeersels'

fnameInfo = strcat(workingDir, '\xOutputs\', fnameBase, '-info.txt'); % infofile name
fInfo = fopen(fnameInfo, 'w'); % create new file with simulation info
tStart = tic; % starting time for timer
CPUtimeAtStart = cputime; % CPU start time

fprintf( fInfo, 'Start: %s\n', datestr(now, 'dd mmm yy @ HH:MM:SS')); % print start time

%% Load data

if dims == 1 % 1D example
    fname = strcat(workingDir, '\xWorkspace-data\', 'GPEpaper2_Example1_1D-workspace.mat');
elseif dims == 2 % 2D example
    fname = strcat(workingDir, '\xWorkspace-data\', 'GPEpaper1_Example1_2D-workspace.mat');
elseif dims == 3 % 3D example
    fname = strcat(workingDir, '\xWorkspace-data\', 'GPEpaper2_Example4_3D-workspace-groundstate.mat'); 
else
    printf('Dimensionality not given!')
end

load(fname)

%% Save the workspace & simulation info

endTime = now;
tElapsed = toc(tStart); % save time elapsed
%dimensions = str2double(fnameBase(20)); % get dimensionality from file name

% save information about final simulation iteration in info file
add_info( fInfo, Outputs, Method, dims, tElapsed, CPUtimeAtStart, cputime); 
fprintf( fInfo, 'End: %s\n', datestr(endTime, 'dd mmm yy @ HH:MM:SS')); % end time
fclose( fInfo );

% save workspace to workspace folder
fnameWspace = strcat('xWorkspace-data/', fnameBase);
%save(strcat(fnameWspace, '-workspace') )
%save(strcat(fnameWspace, '-outputs') , 'Phi', 'Outputs')