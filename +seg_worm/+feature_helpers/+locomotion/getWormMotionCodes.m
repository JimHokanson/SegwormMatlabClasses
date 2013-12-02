function all_events_struct = getWormMotionCodes(midbody_speed, lengths,fps)
%
%   events = seg_worm.feature_helpers.locomotion.getWormMotionCodes(midbody_speed, lengths,fps)
%
%   Output
%   =======================================================================
%   events   - the locomotion events; a struct with event fields:
%
%                  forward  = forward locomotion
%                  paused   = no locomotion (the worm is paused)
%                  backward = backward locomotion
%                  mode     = the locomotion mode:
%
%                             -1 = backward locomotion
%                              0 = no locomotion (the worm is paused)
%                              1 = forward locomotion
%
%
%   Nature Methods Description
%   =======================================================================
%   Motion States. The worm’s forward, backward, and paused motion states
%   attempt to differentiate these event states unambiguously
%   (Supplementary Fig. 4f). Therefore, ambiguous motion has no associated
%   state.
% 
%   The motion states are computed from the worm’s velocity and length
%   (described  in the section on “Morphology”). Missing lengths are
%   linearly interpolated between segmented frames.
%
%   The following filtering criteria were chosen based on human labeling of
%   events within a variety of N2 and mutant videos.
%
%   The worm is defined in a state of forward motion when a period, more
%   than half a second long, is observed wherein:
%           a) the worm travels at least 5% of its mean length over the
%   entire period; and,
%           b) the worm’s speed is at least 5% of its length,
%   per second, in each frame.
%
%   The worm must maintain this speed almost continuously with permissible
%   interruptions of, at most, a quarter second (this permits quick
%   contradictory movements such as head withdrawal, body contractions, and
%   segmentation noise).
%
%   The criteria for backward motion is identical except the worm must be
%   moving backwards (the midbody speed must be negatively signed). 
%
%   The worm is defined in a paused state when a period, more than half a
%   second long, is observed wherein the worm’s forward and backward speed
%   do not exceed 2.5% of its length, per second, in each frame.
%
%   The worm must observe these speed limits almost continuously with
%   permissible interruptions of, at most, a quarter second (once again,
%   this permits quick contradictory movements).



%These are a percentage of the length ...
SPEED_THRESHOLD_PCT   = 0.05;
DISTANCE_THRSHOLD_PCT = 0.05;
PAUSE_THRESHOLD_PCT   = 0.025;

%% Compute the locomotion events.
%--------------------------------------------------------------------------

% Initialize the worm speed and video frames.
totalFrames = length(midbody_speed);

% Compute the distance.
distance = abs(midbody_speed / fps);

% Interpolate the missing lengths.
isNotData = isnan(lengths);
isData    = ~isNotData;
dataI     = find(isData);
interpI   = find(isNotData);
if ~isempty(interpI) && length(dataI) > 1
    lengths(interpI) = interp1(dataI, lengths(dataI), interpI, 'linear');
end

%==========================================================================
wormSpeedThr               = lengths * SPEED_THRESHOLD_PCT;
wormDistanceThr            = lengths * DISTANCE_THRSHOLD_PCT; 
wormEventFramesThr         = 0.5 * fps;
wormEventMinInterFramesThr = 0.25 * fps;

%Forward stuffs
%--------------------------------------------------------------------------
minForwardSpeed            = wormSpeedThr;
minForwardDistance         = wormDistanceThr;

%Backward stuffs
%--------------------------------------------------------------------------
maxBackwardSpeed    = -wormSpeedThr;
minBackwardDistance = wormDistanceThr;

%Paused stuffs
%--------------------------------------------------------------------------
wormPauseThr   = lengths * PAUSE_THRESHOLD_PCT; % 2.5 percent of its length
minPausedSpeed = -wormPauseThr;
maxPausedSpeed = wormPauseThr;

min_speeds   = {minForwardSpeed    []                  minPausedSpeed};
max_speeds   = {[]                 maxBackwardSpeed    maxPausedSpeed};
min_distance = {minForwardDistance minBackwardDistance []};

frames = cell(1,3);

all_events_struct = struct;

FIELD_NAMES = {'forward' 'backward' 'paused'};

motion_mode = NaN(1,totalFrames);

for iType = 1:3
   
    %TODO: This will all be moved into the event class for better
    %encapsulation
    
    %NOTE: Since this is relatively complicated, I might want to create an
    %event finder class which allows specification of options ...
    frames{iType} = seg_worm.feature.event.findEvent( ...
    midbody_speed, ...      1
    min_speeds{iType}, ...  2
    max_speeds{iType}, ...  3
    true, ...               4
    wormEventFramesThr, ... 5
    [], ...                 6
    false, ...              7
    min_distance{iType}, ...8
    [], ...                 9
    true, ...               10
    distance, ...           11
    wormEventMinInterFramesThr); %12
    

    %TODO: Include mode assignment here 


    %TODO: Replace this code with the event class, seems to be ok
    %but I need to clarify the output structure and what forms it takes ...
    
    % Compute the statistics.
    [event_stats, stats] = seg_worm.events.events2stats(frames{iType}, fps, distance, 'distance', 'interDistance');

    % Organize the output ...
    cur_field_name = FIELD_NAMES{iType};
    all_events_struct.(cur_field_name).frames = event_stats;
    if isempty(stats)
        all_events_struct.(cur_field_name).frequency = [];
        all_events_struct.(cur_field_name).ratio     = [];
    else
        all_events_struct.(cur_field_name).frequency = stats.frequency;
        all_events_struct.(cur_field_name).ratio     = stats.ratio;
    end

    wtf = seg_worm.feature.event(frames{iType},fps,distance,'distance','interDistance');    
    wtf2 = wtf.getFeatureStruct;
    keyboard
    
end
all_events_struct.mode = motion_mode;


%% Compute the motion mode. - moving this to be in the above loop
%--------------------------------------------------------------------------
% % % % % Set forward = 1, backward = -1, paused = 0, and unknown = NaN.
% % % % all_events_struct.mode  = nan(1,totalFrames);
% % % % frame_values = [1 -1 0];
% % % % for iType = 1:3
% % % %    %Didn't look at this function, might make static from event ...
% % % %    % Translate the events to logical arrays.
% % % %    mask = seg_worm.events.events2array(frames{iType},  totalFrames);
% % % %    all_events_struct.mode(mask) = frame_values(iType); 
% % % % end


end