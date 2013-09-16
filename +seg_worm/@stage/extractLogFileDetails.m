function [mediaTimes,locations] = extractLogFileDetails(obj)


%1)Real Time        - 21/04/2010 17:19:21.093	
%2)Media Time       - 00:00.0
%3)Location Type    - STAGE
%4)Centroid/Stage/Speed X (microns[/second]) - 11322.7414953351	
%5)Centroid/Stage/Speed Y (microns[/second]) - 16035.0616507836
%
%   NONE OF THESE HAD EXAMPLES ... :/
%   
%6)MER Min X (microns)	
%7)MER Min Y (microns)	
%8)MER Width (microns)	
%9) MER Height (microns)


%FORMAT:
%---------------------------------------
%1) a header line describing each column

logFile = obj.log_file_path;

% Open the log file.
fid = fopen(logFile);
if fid < 0
    error('findStageMovement:BadLogFile', ['Cannot open ' logFile]);
end

% Read the media times and locations from the log file.
mediaTimes = [];
locations = [];
i = 1;
lineNum = 1;
while ~feof(fid)
    
    % Read a line from the log file.
    line = fgetl(fid);
    [values,count] = sscanf(line, ['%*d/%*d/%*d %*d:%*d:%*f,' ...
        '%d:%d:%f,' '%*[^,],' '%f,' '%f,' '%*[^\n]s\n']);
    
    % Extract the media time and location from the line.
    if count == 5
        mediaTime = ((values(1) * 60) + values(2)) * 60 + values(3);
        location = [values(4) values(5)];
        
        % Was the location duplicated?
        if i > 1 && isequal(locations(i - 1,:), location)
            
            % Duplicate locations, prior to recording, are discarded.
            if mediaTime > 0
                error('findStageMovement:DuplicateLocation', ...
                    ['On line ' num2str(lineNum) ', at ' ...
                    num2str(values(1)) ':' ...
                    num2str(values(2), '%02.0f') ':' ...
                    num2str(values(3), '%02.3f') ', location (' ...
                    num2str(location(1)) ', ' num2str(location(2)) ...
                    ') is duplicated.']);
            end
            
        % Advance.
        else
            mediaTimes(i) = mediaTime;
            locations(i,:) = location;
            i = i + 1;
        end
    end
    
    % Advance.
    lineNum = lineNum + 1;
end

% Clean up.
fclose(fid);

end