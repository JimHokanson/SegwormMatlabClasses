classdef event < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.feature.event
    
    properties
       starts
       ends
       times
       inter_times
       inter_names
    end
    
    methods
        function obj = event()
           %TODO: This will be from event2stats
        end
        function getStruct(obj)
            
        end
    end
    
    methods (Static)
        function frames = findEvent(varargin)
           %TODO: Copy over findEvent 
        end
    end
    
end

%{
coilFrames = coilEventStats;
coilFrequency = [];
coilTimeRatio = [];
if ~isempty(coiledStats)
    coilFrequency = coiledStats.frequency;
    coilTimeRatio = coiledStats.ratio.time;
end
posture.coils = struct( ...
    'frames',       coilFrames, ...
    'frequency',    coilFrequency, ...
    'timeRatio',    coilTimeRatio);
%}