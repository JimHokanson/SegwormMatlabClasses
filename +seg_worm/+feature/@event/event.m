classdef event < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.feature.event
    %
    %   See Also:
    %   seg_worm.events.events2stats
    %
    %
    %   General Notes:
    %   -------------------------------------------------------------------
    %   - This is still very much a work in progress. The events are quite
    %   complicated and messy
    %   - An event can optionally contain
    
    %Known Uses:
    %----------------------------------------------------------------------
    %posture.coils - seg_worm.feature_helpers.posture.getCoils
    %
    %locomotion.turns.omegas
    %locomotion.turns.upsilons
    %
    %
    %Uses findEvent ...
    %----------------------------------------------------------------------
    %  in seg_worm.feature_helpers.locomotion.getWormMotionCodes
    %locomotion.motion.forward  
    %locomotion.motion.backward
    %locomotion.motion.paused
    
    %.frames    - event_stats (from event2stats)
    %.frequency -
    
    %Final outputs
    %{
    
    .frames
        .start
        .end
        .time
        .interTime
        .(data_sum_name) - if specified
        .(inter_data_sum_name) - if specified
    .frequency
    %
    %}
    
    properties
        fps
        n_video_frames
        
        %INPUTS
        %------------------------------------------------------------------
        start_Is %[1 n_events]
        end_Is   %[1 n_events]
        data_sum_name %[1 n_events]
        inter_data_sum_name %[1 n_events], last value is NaN
        
        %Outputs - see events2stats
        %------------------------------------------------------------------
        event_durations %[1 n_events]
        inter_event_durations %[1 n_events], last value is NaN
        
        %These two properties are missing if the input names
        %are empty
        data_sum_values
        inter_data_sum_values
        
        total_time
        frequency
        time_ratio
        data_ratio %[1 1] might not exist if .data_sum_name is not specified
        
    end
    
    properties (Dependent)
        n_events
        n_events_for_stats
    end
    
    methods
        function value = get.n_events(obj)
            value = length(obj.start_Is);
        end
        function value = get.n_events_for_stats(obj)
            % Compute the number of events, excluding the partially recorded ones.
            value = obj.n_events;
            if value > 1
                if obj.start_Is(1) == 1
                    value = value - 1;
                end
                if obj.end_Is(end) == obj.n_events
                    value = value - 1;
                end
            end
        end
    end
    
    methods (Static)
       %This is temporary, I'll probably create an event finder class ...
       %
       %    Used by:
       %    locomotion.motion.forward  
       %    locomotion.motion.backward
       %    locomotion.motion.paused
       %
       %    in seg_worm.feature_helpers.locomotion.getWormMotionCodes
       frames = findEvent(data, minThr, maxThr, varargin); 
       function s = getNullStruct(fps,data_sum_name,inter_data_sum_name)
           %
           %
           %    s = seg_worm.feature.event.getNullStruct(fps,data_sum_name,inter_data_sum_name)
           %
           
          event_ss = seg_worm.feature.event_ss([],[]); 
          obj = seg_worm.feature.event(event_ss,fps,[],data_sum_name,inter_data_sum_name);
          s = obj.getFeatureStruct();
       end
    end
    
    methods
        function obj = event(event_ss,fps,data,data_sum_name,inter_data_sum_name)
            %
            %   obj = seg_worm.feature.event(event_ss,fps,data,data_sum_name,inter_data_sum_name)
            %
            %
            %   Inputs
            %   ===========================================================
            %   event_ss : seg_worm.feature.event_ss
            %          .start_Is - frame numbers in which events start
            %          .end_Is   - frame numbers in which events end
            %
            %   fps      : (scalar) frames per second
            %   data     : This data is used for computations, it is
            %               either:
            %             1) distance
            %
            %               From: worm.locomotion.velocity.midbody.speed
            %               distance = abs(speed / fps);
            %
            %               Known users:
            %               seg_worm.feature_helpers.posture.getCoils
            %               seg_worm.feature_helpers.locomotion.getWormMotionCodes
            %
            %             2) ????
            %
            %   data_sum_name : (char) When retrieving the final structure
            %         this is the name given to the field that contains the
            %         sum of the input data during the event
            %   inter_data_sum_name : (char) "          " sum of the input
            %         data between events
            %
            %Some of this code is based on event2stats
            
            obj.fps            = fps;
            obj.n_video_frames = length(data);
            
            if isobject(event_ss)
                obj.start_Is  = event_ss.start_Is;
                obj.end_Is    = event_ss.end_Is;
            else
                obj.start_Is  = [event_ss.start];
                obj.end_Is    = [event_ss.end];                
            end
            obj.data_sum_name       = data_sum_name;
            obj.inter_data_sum_name = inter_data_sum_name;
            
            %Now populate the outputs ...
            %--------------------------------------------------------------
            if obj.n_events == 0
                return
            end

            %---------------------------
            obj.event_durations       = (obj.end_Is - obj.start_Is + 1)./obj.fps;
            obj.inter_event_durations = [obj.start_Is(2:end) - obj.end_Is(1:end-1) - 1 NaN]./fps;
            
            %---------------------------
            if ~isempty(obj.data_sum_name)
                temp = zeros(1,obj.n_events);
                for iEvent = 1:obj.n_events
                    temp(iEvent) = nansum(data(obj.start_Is(iEvent):obj.end_Is(iEvent)));
                end
                obj.data_sum_values = temp;
            end
            
            %---------------------------
            if ~isempty(inter_data_sum_name)
                temp = NaN(1,obj.n_events);
                for iEvent = 1:(obj.n_events-1)
                    start_frame = obj.end_Is(iEvent)+1;
                    end_frame   = obj.start_Is(iEvent+1)-1;
                    temp(iEvent) = nansum(data(start_frame:end_frame));
                end
                obj.inter_data_sum_values = temp;
            end
            
            %----------------------------
            obj.total_time = obj.n_video_frames/obj.fps;
            obj.frequency  = obj.n_events_for_stats/obj.total_time;
            
            obj.time_ratio = nansum(obj.event_durations) / obj.total_time;
            if ~isempty(obj.data_sum_name)
               obj.data_ratio = nansum(obj.data_sum_values)/nansum(data); 
            end
        end
        function s = getFeatureStruct(obj)
            %   
            %   This function returns the structure that matches the form
            %   seen in the feature files
            %
            %   JAH TODO: Describe format of structure ...
            %
            %   
            
            %This bit of code is meant to replace all of the little
            %extra stuff that was previously being done after getting
            %the event frames and converting them to stats
            
            s = struct;
            
            if obj.n_events == 0
               s = struct('frames',[],'frequency',[],'timeRatio',[]);
               return
            end
            
            %--------------------------------------------------------------
            f = struct(...
                'start',        num2cell(obj.start_Is),...
                'end',          num2cell(obj.end_Is), ...
                'time',         num2cell(obj.event_durations),...
                'interTime',    num2cell(obj.inter_event_durations));
            
            if ~isempty(obj.data_sum_name)
               temp = num2cell(obj.data_sum_values);
               [f.(obj.data_sum_name)] = deal(temp{:});
            end
            
            if ~isempty(obj.inter_data_sum_name)
               temp = num2cell(obj.inter_data_sum_values);
               [f.(obj.inter_data_sum_name)] = deal(temp{:});
            end
            
            %This is correct for coiled events, not sure about others ...
            %--------------------------------------------------------------
            s.frames    = f;
            s.frequency = obj.frequency;
            
            
            %??? - why the difference, how to know ????
            %------------------------------------------------
            %ratio struct is present if worm can travel during event
            %
            %  - this might correspond to data_sum being defined 
            %
            %- for motion codes - data and interdata
            %ratio.time
            %ratio.distance
            %
            %- for coils - just interdata
            %timeRatio - no ratio field
            
            %?? Do we also need a check for inter_data as well?
            if isempty(obj.data_sum_name)
                s.timeRatio = obj.time_ratio;
            else
                s.ratio.time     = obj.time_ratio;
                s.ratio.distance = obj.data_ratio;
            end
        end
        function new_obj = mergeEvents(obj1,obj2)
            
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