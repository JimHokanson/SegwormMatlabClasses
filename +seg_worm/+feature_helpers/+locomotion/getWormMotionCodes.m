function events = getWormMotionCodes(lengths,fps)
%
%   seg_worm.feature_helpers.locomotion.getWormMotionCodes(lengths)
%

%JAH: Code unfinished ...

%% Compute the locomotion events.
%--------------------------------------------------------------------------
events = [];

% Initialize the worm speed and video frames.
speed = velocity.midbody.speed;
totalFrames = length(speed);

% Compute the distance.
distance = abs(speed / fps);

% Interpolate the missing lengths.
isNotData = isnan(lengths);
isData = ~isNotData;
dataI = find(isData);
interpI = find(isNotData);
if ~isempty(interpI) && length(dataI) > 1
    lengths(interpI) = interp1(dataI, lengths(dataI), interpI, 'linear');
end



%% Find the forward motion.
%--------------------------------------------------------------------------
wormSpeedThr               = lengths * 0.05; % 5 percent of its length
wormDistanceThr            = lengths * 0.05; % 5 percent of its length
wormEventFramesThr         = 0.5 * fps;
wormEventMinInterFramesThr = 0.25 * fps;
minForwardSpeed            = wormSpeedThr;
minForwardDistance         = wormDistanceThr;
forwardFrames = findEvent(speed, minForwardSpeed, [], true, ...
    wormEventFramesThr, [], false, ...
    minForwardDistance, [], true, distance, wormEventMinInterFramesThr);

% Compute the forward statistics.
[forwardEventStats, forwardStats] = events2stats(forwardFrames, fps, distance, 'distance', 'interDistance');

% Organize the forward motion.
events.forward.frames = forwardEventStats;
events.forward.frequency = [];
events.forward.ratio = [];
if ~isempty(forwardStats)
    events.forward.frequency = forwardStats.frequency;
    events.forward.ratio = forwardStats.ratio;
end



%% Find the backward motion.
maxBackwardSpeed    = -wormSpeedThr;
minBackwardDistance = wormDistanceThr;
backwardFrames = findEvent(speed, [], maxBackwardSpeed, true, ...
    wormEventFramesThr, [], false, ...
    minBackwardDistance, [], true, distance, wormEventMinInterFramesThr);

% Compute the backward statistics.
[backwardEventStats, backwardStats] = events2stats(backwardFrames, fps, ...
    distance, 'distance', 'interDistance');

% Organize the backward motion.
events.backward.frames    = backwardEventStats;
events.backward.frequency = [];
events.backward.ratio     = [];
if ~isempty(backwardStats)
    events.backward.frequency = backwardStats.frequency;
    events.backward.ratio = backwardStats.ratio;
end



%% Find the paused motion.
wormPauseThr   = lengths * 0.025; % 2.5 percent of its length
minPausedSpeed = -wormPauseThr;
maxPausedSpeed = wormPauseThr;
pausedFrames = findEvent(speed, minPausedSpeed, maxPausedSpeed, true, ...
    wormEventFramesThr, [], false, ...
    [], [], true, distance, wormEventMinInterFramesThr);

% Compute the paused statistics.
[pausedEventStats, pausedStats] = events2stats(pausedFrames, fps, distance, 'distance', 'interDistance');

% Organize the paused motion.
events.paused.frames    = pausedEventStats;
events.paused.frequency = [];
events.paused.ratio     = [];
if ~isempty(pausedStats)
    events.paused.frequency = pausedStats.frequency;
    events.paused.ratio = pausedStats.ratio;
end



%% Compute the motion mode.

% Translate the events to logical arrays.
isForwardFrame  = events2array(forwardFrames, totalFrames);
isBackwardFrame = events2array(backwardFrames, totalFrames);
isPausedFrame   = events2array(pausedFrames, totalFrames);

% Set forward = 1, backward = -1, paused = 0, and unknown = NaN.
events.mode = nan(1,totalFrames);
events.mode(isForwardFrame) = 1;
events.mode(isBackwardFrame) = -1;
events.mode(isPausedFrame) = 0;

end