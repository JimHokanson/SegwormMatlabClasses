classdef file_name_info < handle
    %
    %   Class:
    %   seg_worm.video.file_name_info
    %
    
    properties
       file_name 
       strain
       food
       side
       timestamp
    end
    
    methods
        function obj = file_name_info(avi_file_path)
        %
        %
        %   obj = seg_worm.video.file_name_info(avi_file_path)
        
        
        %This might be the wrong file name, side is not in the video file, perhaps this
        %is added later ?????
            
           [~,obj.file_name] = fileparts(avi_file_path);
            
pattern = [...
    '(?<type>\S+)\s*' ...
    '(?<strain>\(\w+\)\w*)?\s+' ...
    '(?<food>on\s+food)?\s*' ...
    '(?<side>[LR])_' ...
    '(?<year>\d\d\d\d)_' ...
    '(?<month>\d\d)_' ...
    '(?<day>\d\d)__' ...
    '(?<hour>\d\d)_' ...
    '(?<minute>\d\d)_' ...
    '(?<second>\d\d).*'];
matches = regexpi(obj.file_name, pattern, 'names');
if isempty(matches)
    error('parseWormFilename:NoMatches', 'The filename does not parse');
end

% Parse the strain.
smatches(1) = struct('strain', [], 'chromosome', []);
if ~isempty(matches(1).strain)
    spattern = '\((?<strain>\w+)\)(?<chromosome>\w+)?';
    smatches = regexpi(matches(1).strain, spattern, 'names');
end
    
% Construct the date.
date = datenum(str2num(matches(1).year), ...
    str2num(matches(1).month), ...
    str2num(matches(1).day), ...
    str2num(matches(1).hour), ...
    str2num(matches(1).minute), ...
    str2num(matches(1).second));
    
% Fill the info.
info = struct('type', matches(1).type, ...
    'strain',         smatches(1).strain, ...
    'chromosome',     smatches(1).chromosome, ...
    'food',          ~isempty(matches(1).food), ...
    'side',           matches(1).side, ...
    'timestamp',      date);
            
        end
    end
    
end

