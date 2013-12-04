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
%
%   Questions:
%   - based on the way these are calculated it seems like a frame
%   could have no event type (i.e. not be going forward, backward, or
%   paused)


%These are a percentage of the length ...
SPEED_THRESHOLD_PCT   = 0.05;
DISTANCE_THRSHOLD_PCT = 0.05;
PAUSE_THRESHOLD_PCT   = 0.025;

%% Compute the locomotion events.
%--------------------------------------------------------------------------

% Initialize the worm speed and video frames.
n_frames = length(midbody_speed);

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

%--------------------------------------------------------------------------

all_events_struct = struct;

FIELD_NAMES  = {'forward' 'backward' 'paused'};
FRAME_VALUES = [1 -1 0];
motion_mode = NaN(1,n_frames);

for iType = 1:3
   
    %Determine when the event type occurred
    %----------------------------------------------------------------------
    ef = seg_worm.feature.event_finder;
    
    ef.include_at_thr = true;
    ef.min_frames_thr = wormEventFramesThr;
    ef.min_sum_thr    = min_distance{iType};
    ef.include_at_sum_thr   = true;
    ef.data_for_sum_thr     = distance;
    ef.min_inter_frames_thr = wormEventMinInterFramesThr;
    
    frames_temp = ef.getEvents(midbody_speed,min_speeds{iType},max_speeds{iType});
    
    %Assign event type to relevant frames
    %----------------------------------------------------------------------
    mask = frames_temp.getEventMask(n_frames);
    motion_mode(mask) = FRAME_VALUES(iType);

    %Take the start and stop indices and convert them to the structure
    %used in the feature files ...
    %----------------------------------------------------------------------
    cur_field_name = FIELD_NAMES{iType};

    temp = seg_worm.feature.event(frames_temp,fps,distance,'distance','interDistance');    
    all_events_struct.(cur_field_name) = temp.getFeatureStruct;
    
end
all_events_struct.mode = motion_mode;


end