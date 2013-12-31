function getWormMotionCodes(obj,midbody_speed, skeleton_lengths,fps)
%
%   events = seg_worm.feature_helpers.locomotion.getWormMotionCodes(midbody_speed, lengths,fps)
%
%   Inputs
%   =======================================================================
%   midbody_speed    : [1 x n_frames] from locomotion.velocity.midbody.speed
%   skeleton_lengths : [1 x n_frames]
%   fps              : (scalar) frames per second
%
%
%   Output
%   =======================================================================
%   all_events_struct : the locomotion events; a struct with event fields:
%
%              forward  - (event) forward locomotion
%              paused   - (event) no locomotion (the worm is paused)
%              backward - (event) backward locomotion
%              mode     = [1 x n_frames] the locomotion mode:
% 
%                         -1 = backward locomotion
%                          0 = no locomotion (the worm is paused)
%                          1 = forward locomotion
%
%
%   Old Name: 
%   - part of wormVelocity.m
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

%These are times (s)
EVENT_FRAMES_THRESHOLD = 0.5; %Half a second
EVENT_MIN_INTER_FRAMES_THRESHOLD = 0.25;

DATA_SUM_NAME = 'distance';
INTER_DATA_SUM_NAME = 'interDistance';


%% Compute the locomotion events.
%--------------------------------------------------------------------------

% Initialize the worm speed and video frames.
n_frames = length(midbody_speed);

% Compute the distance.
distance = abs(midbody_speed / fps);

% Interpolate the missing lengths.
isNotData = isnan(skeleton_lengths);
isData    = ~isNotData;
dataI     = find(isData);
interpI   = find(isNotData);
if ~isempty(interpI) && length(dataI) > 1
    skeleton_lengths(interpI) = interp1(dataI, skeleton_lengths(dataI), interpI, 'linear');
end

%==========================================================================
worm_speed_threshold       = skeleton_lengths * SPEED_THRESHOLD_PCT;
worm_distance_threshold    = skeleton_lengths * DISTANCE_THRSHOLD_PCT; 

%Forward stuffs
%--------------------------------------------------------------------------
min_forward_speed    = worm_speed_threshold;
min_forward_distance = worm_distance_threshold;

%Backward stuffs
%--------------------------------------------------------------------------
max_backward_speed    = -worm_speed_threshold;
min_backward_distance = worm_distance_threshold;

%Paused stuffs
%--------------------------------------------------------------------------
worm_pause_threshold = skeleton_lengths * PAUSE_THRESHOLD_PCT; % 2.5 percent of its length
min_paused_speed     = -worm_pause_threshold;
max_paused_speed     = worm_pause_threshold;

min_speeds   = {min_forward_speed       []                      min_paused_speed};
max_speeds   = {[]                      max_backward_speed      max_paused_speed};
min_distance = {min_forward_distance    min_backward_distance   []};

%--------------------------------------------------------------------------
worm_event_frames_threshold          = fps * EVENT_FRAMES_THRESHOLD;
worm_event_min_interframes_threshold = fps * EVENT_MIN_INTER_FRAMES_THRESHOLD;

all_events_struct = struct;

FIELD_NAMES  = {'forward' 'backward' 'paused'};
FRAME_VALUES = [1 -1 0];
motion_mode = NaN(1,n_frames);

for iType = 1:3
   
    %Determine when the event type occurred
    %----------------------------------------------------------------------
    ef = seg_worm.feature.event_finder;
    
    ef.include_at_thr = true;
    ef.min_frames_thr = worm_event_frames_threshold;
    ef.min_sum_thr    = min_distance{iType};
    ef.include_at_sum_thr   = true;
    ef.data_for_sum_thr     = distance;
    ef.min_inter_frames_thr = worm_event_min_interframes_threshold;
    
    %seg_worm.feature.event_finder.getEvents
    frames_temp = ef.getEvents(midbody_speed,min_speeds{iType},max_speeds{iType});
    %frames_temp - class - seg_worm.feature.event_ss
    
    %Assign event type to relevant frames
    %----------------------------------------------------------------------
    mask = frames_temp.getEventMask(n_frames);
    motion_mode(mask) = FRAME_VALUES(iType);

    %Take the start and stop indices and convert them to the structure
    %used in the feature files ...
    %----------------------------------------------------------------------
    cur_field_name = FIELD_NAMES{iType};

    temp = seg_worm.feature.event(frames_temp,fps,distance,DATA_SUM_NAME,INTER_DATA_SUM_NAME);    
    all_events_struct.(cur_field_name) = temp.getFeatureStruct;
    
end

all_events_struct.mode = motion_mode;

obj.motion = all_events_struct;

end