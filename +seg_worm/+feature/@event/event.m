classdef event < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.feature.event
    
    properties
        fps
        n_samples
        starts %[1 n_samples]
        ends   %[1 n_samples]
        times  %????
        inter_times    %[1 n_samples]
        
        inter_distance %
        %posture.coils
        %locomotion.turns.omegas
        %locomotion.turns.upsilons
        
        is_ventral %
        %locomotion.turns.omegas
        %locomotion.turns.upsilons
        
        inter_names  %[1 n_samples]
    end
    
    properties (Dependent)
        n_full_events
    end
    
    methods
        function value = get.n_full_events(obj)
            % Compute the number of events, excluding the partially recorded ones.
            value = length(obj.starts);
            if value > 1
                if obj.starts == 1
                    value = value - 1;
                end
                if obj.ends(end) == obj.n_samples
                    value = value - 1;
                end
            end
        end
    end
    
    methods
        function obj = event()
            %TODO: This will be from event2stats
        end
        function getStruct(obj)
            %This will be to get the structure for saving ...
            
            s = struct;
            f = struct;
            
            %TODO: What does this look like ???
            
            keyboard
            
            s.frames = struct(...
                'start',num2cell(obj.starts),...
                'end',num2cell(obj.ends));
            
            %FORMAT 1 (coils)
            %--------------------------------
            %.
            %.frequency
            %.timeRatio
            
            
            
            
            % Compute the event statistics.
            %--------------------------------------------------------------
            totalTime = length(data) / fps;
            frequency = numEvents / totalTime;
            timeRatio = nansum([eventStats.time]) / totalTime;
            if isName
                dataRatio = nansum([eventStats.(name)]) / nansum(data);
                ratios = struct( ...
                    'time',     timeRatio, ...
                    name,       dataRatio);
            else
                ratios = struct('time', timeRatio);
            end
            summaryStats = struct( ...
                'frequency', frequency, ...
                'ratio',     ratios);
            
            
            % %Reorganize everything for the feature file.
            % omegaFrequency = [];
            % omegaTimeRatio = [];
            % if ~isempty(omegaStats)
            %     omegaFrequency = omegaStats.frequency;
            %     omegaTimeRatio = omegaStats.ratio.time;
            % end
            % omegas = struct( ...
            %     'frames', omegaFrames, ...
            %     'frequency', omegaFrequency, ...
            %     'timeRatio', omegaTimeRatio);
            
            
            
        end
        function new_obj = mergeEvents(obj1,obj2)
            
        end
    end
    
    methods (Static)
        %This is in another file, I want this to lead
        %to
        frames = findEvent(data, minThr, maxThr, varargin)
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