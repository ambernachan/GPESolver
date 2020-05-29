
classdef Info  < handle
    properties
        name
        startTime
        startTimeString
        startCpuTime
        timerStartValue

        outputdir
        subdir
        fulldir
        filename
        path
        fileID
    end
    methods
        % 
        function obj = Info(name)
            if nargin ~= 1
                error('Invalid number of arguments given %s, expected 1', nargin);
            end
            if ~(isstring(name) || ischar(name))
                error('Name should be a string or character array');
            end
            obj.name = string(name);
            obj.startTime = now;
            obj.startTimeString = datestr(obj.startTime, 'yyyy-dd-mm_HH-MM-SS');
            obj.startCpuTime = cputime;
            obj.timerStartValue = tic

            % Set dirs and filenames
            obj.outputdir = 'xOutputs';
            obj.subdir = name;
            obj.fulldir = string(pwd) + '/' + obj.outputdir + '/' + obj.subdir;
            obj.filename = string('info_') + obj.startTimeString + '.txt';
            obj.path = obj.fulldir + '/' + obj.filename;
            
            % Create directory
            [mkdirStatus, mkdirMsg] = mkdir(char(obj.fulldir));
            if ~mkdirStatus
                error('Unable to create directory: %s. Error: %s', obj.fulldir, mkdirMsg);
            end
            
            % Open info file
            [obj.fileID, errmsg] = fopen(obj.path, 'w');
            if obj.fileID == -1
                error('Unable to open file: %s', errmsg);
            end 

            % Initially print the start time
            fprintf(obj.fileID, 'Start: %s\n', obj.startTimeString); 
        end
        
        function add_info(obj)
            if obj.fileID == -1
                error('File not open, did you already finish this Info object?');
            end
        end
        
        function finish(obj)
            if obj.fileID == -1
                error('File not open, did you already finish this Info object?');
            end
            
            fclose(obj.fileID);
            obj.fileID = -1;
        end
    end
end
