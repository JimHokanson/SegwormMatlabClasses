classdef event_finder < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.feature.event_finder
    %   
    %
    %   TODO: Return class of event_ss
    %
    %   See Also:
    %   seg_worm.feature.event_ss
    %   seg_worm.feature.event
    %
    %   Usage:
    %   1) Call constructor to get object for setting options
    %   2) Call getEvents() to get data
    %
    
    properties
       %OPTIONS
       include_at_thr = false %(isAtThr)
       
       %---------------------------------------
       %Sample based
       min_frames_thr = [] %(minFramesThr)
       max_frames_thr = []
       include_at_frames_thr = false
       
       %---------------------------------------
       %Data based
       min_sum_thr = []
       max_sum_thr = []
       include_at_sum_thr = false
       
       %----------------------------------------
       data_for_sum_thr = [] %(sumData) the data
       %for thresholding based on the sum, if empty the event
       %data is used
       
       %----------------------------------------
       %Sample based
       min_inter_frames_thr = []
       max_inter_frames_thr = [] %??? When would anyone want
       %to join events if the time between them is too long???
       include_at_inter_frames_thr = false
       
       %------------------------------------------
       %JAH: This won't work, old code didn't support it
       %Data based
       min_inter_sum_thr = []
       max_inter_sum_thr = []
       include_at_inter_sum_thr = false
       
       %????
       %- some of these parameters seem like they could be contradictory
       %    TODO: I'll need to keep an eye out for this in the code ...
       
    end
    
    methods
       %MAIN FUNCTION:
       %seg_worm.feature.event_finder.getEvents
    end
    
end

