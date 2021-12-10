function newVariable = importVariable(workspaceName, variableName)
%IMPORTFILE(workspaceName, variableName)
%  Imports data from the specified file
%  workspaceName:  file to read
%  variableName: variable to import from workspace (optional)

if nargin < 2
    sprintf('Give a variable name to import');
    return
end

% Import the file
workspaceData = load('-mat', workspaceName);

% Create new variables in the base workspace from those fields.
vars = fieldnames(workspaceData);
newVariable = [];

for i = 1:length(vars)
    if(strcmp(vars{i}, variableName))
        %assignin('base', vars{i}, workspaceData.(vars{i}));
        newVariable = workspaceData.(vars{i});
        break
    end
end

