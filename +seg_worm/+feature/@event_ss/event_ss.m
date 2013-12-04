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
    
    properties (Dependent)
        n_events
    end
    
    methods
        function value = get.n_events(obj)
            value = length(obj.start_Is);
        end
    end
    
    methods
        function obj = event_ss(start_Is,end_Is)
            %
            %    obj = seg_worm.feature.event_ss(start_Is,end_Is)
            
            if size(start_Is,1) > 1
                start_Is = start_Is';
            end
            
            if size(end_Is,1) > 1
                end_Is = end_Is';
            end
            
            obj.start_Is = start_Is;
            obj.end_Is   = end_Is;
        end
        %seg_worm.events.events2stats - move here
        %fromStruct - from the old struct version ...
        function mask = getEventMask(obj,n_frames)
            %
            %
            %     Old Name: events2array.m
            
            mask = false(n_frames,1);
            start_Is = obj.start_Is;
            end_Is   = obj.end_Is;
            for iFrame = 1:obj.n_events
                mask(start_Is(iFrame):end_Is(iFrame)) = true;
            end
            
            
        end
    end
    
end

