
classdef Info  < handle
    properties
        name
        dimensions

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
        function obj = Info(name, dimensions)
            if nargin ~= 2
                error('Invalid number of arguments given %d, expected 2', nargin);
            end
            if ~(isstring(name) || ischar(name))
                error('Name should be a string or character array');
            end

            obj.name = string(name);
            obj.dimensions = dimensions;

            obj.separatorStr = '-------------------------------------------';

            obj.creationTime = now;
            obj.creationTimeString = [datestr(obj.creationTime, 'yyyy-mm-dd') '@' datestr(obj.creationTime, 'HH.MM.SS') ];
            obj.creationCpuTime = cputime;
            obj.timerCreationValue = tic;

            % Set dirs and filenames
            obj.outputdir = 'xOutputs';
            obj.subdir = name; % Name of simulation (folder)
            obj.fulldir = string(pwd) + '/' + obj.outputdir + '/' + obj.subdir + '/' + obj.creationTimeString;
            obj.filenameInfo = string('INFO') + '.txt';
            obj.pathInfo = obj.fulldir + '/' + obj.filenameInfo;

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

        function start_info(obj)
            obj.startTime = now;
            obj.startCpuTime = cputime;
            obj.timerStartValue = tic;
        end

        function add_simulation_info(obj, Geometry)
            fprintf(obj.fileID, 'Nx:\t%d\n', Geometry.Nx);
            if obj.dimensions > 1
               fprintf(obj.fileID, 'Ny:\t%d\n', Geometry.Ny);
            end
            if obj.dimensions > 2
               fprintf(obj.fileID, 'Nz:\t%d\n', Geometry.Nz);
            end
        end
        
        % Add basic information about the simulation to an info file
        function add_result_info(obj, Method, Outputs)
            if obj.fileID == -1
                error('File not open, did you already finish this Info object?');
            end
            add_info(obj.fileID, obj.dimensions, Method, Outputs, toc(obj.timerStartValue), obj.startCpuTime, cputime, obj.separatorStr)
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
        function save_workspace(obj, name)
            % save workspace to workspace folder with name = 'groundstate' or 'dynamics'
            workspacePath = obj.fulldir + '/workspace_' + name + '.mat';
            save(char(workspacePath));
        end

        % Save a specified figure (num) to the outputs directory of this simulation
        function save_figure(obj, fignum, state, title)
            figurePath = obj.fulldir + '/' + state + '_fig_' + title + '.fig';
            savefig(fignum, char(figurePath))
        end

        function add_info_separator(obj)
            fprintf(obj.fileID, '%s\n', obj.separatorStr);
        end
    end
end
