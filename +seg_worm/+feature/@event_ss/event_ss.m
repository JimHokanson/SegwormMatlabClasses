classdef event_ss < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.feature.event_ss
    %   
    %   I was going to leave this class as just a Matlab structure but
    %   there is at least one function that would be better as a method
    %   of this class
    %
    %   See Also:
    %   seg_worm.feature.event_finder
    %   seg_worm.feature.event
    
    properties
       start_Is %[1 x n_events]
       end_Is   %[1 x n_events]
    end
    
    methods
        function obj = event_ss(start_Is,end_Is)
           %
           %    obj = seg_worm.feature.event_ss(start_Is,end_Is)
           
           obj.start_Is = start_Is;
           obj.end_Is   = end_Is;
        end
       %seg_worm.events.events2stats - move here
       %fromStruct - from the old struct version ...
    end
    
end

