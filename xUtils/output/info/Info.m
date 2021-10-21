
classdef Info  < handle
    properties
        name
        params

        separatorStr

        creationTime
        creationTimeString
        creationCpuTime
        timerCreationValue

        startTime
        endTime
        startCpuTime
        endCpuTime
        timerStartValue
        
        outputdir
        subdir
        fulldir
        filenameInfo
        pathInfo
        fileID
    end
    methods
        % Constructor
        function obj = Info(name, creationTime, params)
            if nargin ~= 3
                error('Invalid number of arguments given %d, expected 3', nargin);
            end
            if ~(isstring(name) || ischar(name))
                error('Name should be a string or character array');
            end
            if isstring(name)
                name = char(name);
            end

            obj.name = name;
            obj.params = params;

            obj.separatorStr = '-------------------------------------------';

            obj.creationTime = creationTime;
            obj.creationTimeString = [datestr(obj.creationTime, 'yyyy-mm-dd') '@' datestr(obj.creationTime, 'HH.MM.SS') ];
            obj.creationCpuTime = cputime;
            obj.timerCreationValue = tic;

            % Set dirs and filenames
            obj.outputdir = '../xOutputs';
            obj.subdir = name; % Name of simulation (folder)
            obj.fulldir = [pwd '/' obj.outputdir '/' obj.subdir '/' obj.creationTimeString];
            if isfield(obj.params, 'S')
                obj.filenameInfo = ['INFO' obj.suffix() '.txt'];
            else
                obj.filenameInfo = 'INFO.txt';
            end
            obj.pathInfo = [obj.fulldir '/' obj.filenameInfo];

            % Create directory
            [mkdirStatus, mkdirMsg] = mkdir(char(obj.fulldir));
            if ~mkdirStatus
                error('Unable to create directory: %s. Error: %s', obj.fulldir, mkdirMsg);
            end
            
            % Open info file
            [obj.fileID, errmsg] = fopen(obj.pathInfo, 'w');
            if obj.fileID == -1
                error('Unable to open file: %s', errmsg);
            end 

            % Initially print the start time
            fprintf(obj.fileID, 'Start: %s\n', datestr(obj.creationTime, 'dd mmm yy @ HH:MM:SS')); 
            obj.start_info()
        end
        
        function s = suffix(obj)
            if isfield(obj.params,'S')
                if ~isempty(obj.params.S)
                    s = ['_S=' format_str(obj.params.S)];
                else
                    s = '';
                end
            else
                s = '';
            end
        end

        function start_info(obj)
            obj.startTime = now;
            obj.startCpuTime = cputime;
            obj.timerStartValue = tic;
        end

        function add_simulation_info(obj, Geometry)
            fprintf(obj.fileID, 'Nx:\t%d\t\t-%d<x<%d\n', Geometry.Nx, (Geometry.Lx/2), (Geometry.Lx/2));
            if obj.params.dimensions > 1
               fprintf(obj.fileID, 'Ny:\t%d\t\t-%d<y<%d\n', Geometry.Ny, (Geometry.Ly/2), (Geometry.Ly/2));
            end
            if obj.params.dimensions > 2
               fprintf(obj.fileID, 'Nz:\t%d\t\t-%d<z<%d\n', Geometry.Nz, (Geometry.Lz/2), (Geometry.Lz/2));
            end
        end
        
        % Add basic information about the simulation to an info file
        function add_result_info(obj, Method, Outputs)
            if obj.fileID == -1
                error('File not open, did you already finish this Info object?');
            end
            add_info(obj.fileID, obj.params.dimensions, Method, Outputs, toc(obj.timerStartValue), obj.startCpuTime, cputime, obj.separatorStr)
            fprintf(obj.fileID, 'End simulation stage: %s\n', datestr(now, 'dd mmm yy @ HH:MM:SS')); % end time
        end
        
        % Add a string to the INFO.txt file
        function add_custom_info(obj, format_str, varargin)
            fprintf(obj.fileID, format_str, varargin{:});
        end
      
        % Finish up this object, close files and write final information
        function finish_info(obj)
            if obj.fileID == -1
                error('File not open, did you already finish this Info object?');
            end

            obj.endTime = now;
            obj.endCpuTime = cputime;

            obj.add_info_separator()
            fprintf(obj.fileID, 'Total CPU time:\t%8.2f\n', obj.endCpuTime - obj.creationCpuTime);
            fprintf(obj.fileID, 'Total elapsed time:\t' );
            fprintf(obj.fileID, print_time(toc(obj.timerCreationValue)) );
            obj.add_info_separator()

            fprintf(obj.fileID, 'End: %s\n', datestr(obj.endTime, 'dd mmm yy @ HH:MM:SS'));
            
            fclose(obj.fileID);
            obj.fileID = -1;
        end

        % Save a named workspace snapshot to disk
        function wspath = get_workspace_path(obj, name)
            % save workspace to workspace folder with name = 'groundstate' or 'dynamics'
            wspath = [obj.fulldir '/workspace_' name obj.suffix() '.mat'];
            % save(workspacePath);
        end

        % Save a specified figure (num) to the outputs directory of this simulation
        function save_figure(obj, fignum, state, title, path, extension)
            if nargin < 5 || ~exist('extension','var') || isempty(extension)
                extension = '.fig';
                if nargin < 4 || ~exist('path','var') || isempty(path)
                    path = obj.fulldir;
                end
            end
            if ~exist(path, 'dir')
                % path doesn't exist
                % Either prefix is D:\surfdrive\UNI or
                % C:\AMBER\surfdrive\UNI
                elements = strsplit(path, {'/', '\'});
                for i = 1:numel(elements)
                    if strcmp(elements{i}, 'surfdrive')
                        idx = i;
                        return;
                    end
                end
                newpathsuffix = strjoin(elements{idx:end}, '\');
                prefixes = [{'C:\AMBER\'}, {'D:\'}];
                if strcmp(elements{1}, 'C:')
                    prefix = prefixes{2};
                elseif strcmp(elements{1}, 'D:')
                    prefix = prefixes{1};
                end

                path = [prefix newpathsuffix];
            end
            
            if ~exist('title','var') || isempty(title)
                figurePath = [path '/' state obj.suffix() extension];
            else
                figurePath = [path '/' state '_' title obj.suffix() extension];
            end
            %savefig(fignum, char(figurePath))
            
            if strcmp('.fig',extension)
                savefig(fignum, figurePath);
            else
                figHandle = findobj( 'Type', 'Figure', 'Number', fignum );
                saveas(figHandle, figurePath)
            end
        end

        function add_info_separator(obj)
            fprintf(obj.fileID, '%s\n', obj.separatorStr);
        end
    end
end
