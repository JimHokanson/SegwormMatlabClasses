function events = getWormMotionCodes(midbody_speed, lengths,fps)
%
%   seg_worm.feature_helpers.locomotion.getWormMotionCodes()
%
%
%   %       events   - the locomotion events; a struct with event fields:
%
%                  forward  = forward locomotion
%                  paused   = no locomotion (the worm is paused)
%                  backward = backward locomotion
%                  mode     = the locomotion mode:
%
%                             -1 = backward locomotion
%                              0 = no locomotion (the worm is paused)
%                              1 = forward locomotion

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

events = struct;

FIELD_NAMES = {'forward' 'backward' 'paused'};

for iType = 1:3
   
    %TODO: This will all be moved into the event class for better
    %encapsulation
    
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
    
    % Compute the statistics.
    [event_stats, stats] = seg_worm.events.events2stats(frames{iType}, fps, distance, 'distance', 'interDistance');

    % Organize the output ...
    cur_field_name = FIELD_NAMES{iType};
    events.(cur_field_name).frames = event_stats;
    if isempty(stats)
        events.(cur_field_name).frequency = [];
        events.(cur_field_name).ratio     = [];
    else
        events.(cur_field_name).frequency = stats.frequency;
        events.(cur_field_name).ratio     = stats.ratio;
    end

end

%% Compute the motion mode.
%--------------------------------------------------------------------------
% Set forward = 1, backward = -1, paused = 0, and unknown = NaN.
events.mode  = nan(1,totalFrames);
frame_values = [1 -1 0];
for iType = 1:3
   %Didn't look at this function, might make static from event ...
   % Translate the events to logical arrays.
   mask = seg_worm.events.events2array(frames{iType},  totalFrames);
   events.mode(mask) = frame_values(iType); 
end


end