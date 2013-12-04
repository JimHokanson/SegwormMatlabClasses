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
   

    %JAH TODO: I'm currently working on the getEvents method ...

    ef = seg_worm.feature.event_finder;
    
    ef.include_at_thr = true;
    ef.min_frames_thr = wormEventFramesThr;
    ef.min_sum_thr    = min_distance{iType};
    ef.include_at_sum_thr   = true;
    ef.data_for_sum_thr     = distance;
    ef.min_inter_frames_thr = wormEventMinInterFramesThr;
    
    frames_temp = ef.getEvents(midbody_speed,min_speeds{iType},max_speeds{iType});
    

    %This following bit is old code which the code above is replacing
    %----------------------------------------------------------------------
    %JAH NOTE: Current goal is to verify that frames_temp is equal
    %to frames{iType}
    
    
%       ISATTHR,            4
%       MINFRAMESTHR,       5
%       MAXFRAMESTHR,       6
%       ISATFRAMESTHR,      7
%       MINSUMTHR,          8
%       MAXSUMTHR,          9
%       ISATSUMTHR,         10
%       SUMDATA,            11
%       MININTERFRAMESTHR,  12
%       MAXINTERFRAMESTHR,  13
%       ISATINTERFRAMESTHR, 14
%       MININTERSUMTHR,     15
%       MAXINTERSUMTHR,     16
%       ISATINTERSUMTHR)    17
    
   
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
    
    %End of old code
    %----------------------------------------------------------------------
    
    cur_field_name = FIELD_NAMES{iType};

    temp = seg_worm.feature.event(frames_temp,fps,distance,'distance','interDistance');    
    all_events_struct.(cur_field_name) = temp.getFeatureStruct;
    
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