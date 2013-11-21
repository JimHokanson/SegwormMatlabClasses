function velocity = getWormVelocity(x, y, fps, ventralMode)
%getWormVelocity  Compute the worm velocity (speed & direction) at the
%head-tip/head/midbody/tail/tail-tip
%
%   velocity = seg_worm.feature_helpers.locomotion.getWormVelocity(x, y, fps, ventralMode)
%
%   [VELOCITY EVENTS] = WORMVELOCITY(X, Y, FPS, ventralMode)
%
%   Inputs:
%       x           - the worm skeleton's x-axis coordinates
%       y           - the worm skeleton's y-axis coordinates
%       fps         - the frames/seconds
%       ventralMode - the ventral side mode:
%
%                     0 = the ventral side is unknown
%                     1 = the ventral side is clockwise
%                     2 = the ventral side is anticlockwise
%
%   Outputs:
%       velocity - the worm velocity; each field has subfields for the
%                  "speed" and "direction":
%
%                  headTip = the tip of the head (1/12 the worm at 0.25s)
%                  head    = the head (1/6 the worm at 0.5s)
%                  midbody = the midbody (2/6 the worm at 0.5s)
%                  tail    = the tail (1/6 the worm at 0.5s)
%                  tailTip = the tip of the tail (1/12 the worm at 0.25s)
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

%UNFINISHED CODE ...

if ~exist('ventralMode','var')
    ventralMode = 0;
end

SI = seg_worm.skeleton_indices;
TIP_DIFF  = 0.25;
BODY_DIFF = 0.5;

% Compute the tail-to-head direction.
bodyI     = SI.BODY_INDICES;

diffX     = nanmean(diff(x(bodyI,:), 1, 1), 1);
diffY     = nanmean(diff(y(bodyI,:), 1, 1), 1);
bodyAngle = atan2(diffY, diffX).*(180 / pi);

% Compute the velocity.
%--------------------------------------------------------------------------
FIELD_NAMES       = {'headTip'           'head'          'midbody'       'tail'          'tailTip'};
INDICES_LIST      = {SI.HEAD_TIP_INDICES SI.HEAD_INDICES SI.MID_INDICES SI.TAIL_INDICES SI.TAIL_TIP_INDICES};
TIME_SCALE_VALUES = [TIP_DIFF BODY_DIFF BODY_DIFF BODY_DIFF TIP_DIFF];

for iField = 1:length(FIELD_NAMES)
   cur_field_name = FIELD_NAMES{iField};
   cur_indices    = INDICES_LIST{iField};
   cur_scale      = TIME_SCALE_VALUES(iField);
   velocity.(cur_field_name)  = h__computeVelocity(x, y, bodyAngle, cur_indices, fps, cur_scale, ventralMode);
end

end



%% Compute the velocity.
function velocity = h__computeVelocity(x, y, bodyAngle, pointsI, fps, scale_time, ventralMode)
%
%
%   The velocity is computed not using the nearest values but values
%   that are separated by a sufficient time (scale_time). If the encountered 
%   values are not valid, the width of time is expanded up to a maximum of
%   2*scale_time (or technically, 1 scale time in either direction)
%

n_frames = size(x,2);

%What is scale ????

% The scale must be odd.
%--------------------------------------------------------------
scale_samples_input = scale_time * fps;

scale_low  = floor(scale_samples_input);
scale_high = ceil(scale_samples_input);

if scale_low == scale_high
    if mod(scale_low,2) == 0
        scale_samples_final = scale_samples_input + 1;
    else
        scale_samples_final = scale_low;
    end
elseif mod(scale_high,2) == 0
    scale_samples_final = scale_low;
else
    scale_samples_final = scale_high;
end

speed     = nan(1,n_frames);
direction = nan(1,n_frames);
bodyDirection = nan(1,n_frames);
% Do we have enough coordinates?
if scale_samples_final > n_frames
    velocity.speed     = speed;
    velocity.direction = direction;
    return;
end

% Compute the body part direction.
dX = nanmean(diff(x(pointsI,:), 1, 1), 1);
dY = nanmean(diff(y(pointsI,:), 1, 1), 1);

pointAngle = atan2(dY, dX) * (180 / pi);

% Compute the coordinates.
x_mean = mean(x(pointsI,:), 1);
y_mean = mean(y(pointsI,:), 1);

% Compute the speed using back/front nearest neighbors bounded at twice the scale.
scale_minus_1 = scale_samples_final - 1;
half_scale    = scale_minus_1 / 2;

start_index = 1+half_scale;
end_index   = n_frames-half_scale;

middle_indices     = (start_index:end_index) + half_scale;
%   our start_index frame can only have one valid frame (frame 1)
%   afterwards it is invalid
%
%   To get around this we'll pad with bad values (which we'll never match)
%   then shift the final indices
%
%   e.g.
%   scale = 5
%   start = 3 - this means we grab +/- 2, so we need to have our indices
%       padded by 2 when we grab, or alternatively do a filter (i.e. if statement
%   and removal of indices, which I think is slower)
%
%   NaN NaN 1 2 3
%   1   2   3 4 5  <- temporary indices
%
%   Then at the end shift everything by 2 
%
%   This allows us to remove conditionals ...

matched_left_mask  = false(1,length(middle_indices));
matched_right_mask = false(1,length(middle_indices));

%This tells us whether each value is useable or not for velocity
%Better to do this out of the loop ...
is_good_value_mask = [false(1,half_scale) ~isnan(x_mean) false(1,half_scale)];

%These will be the final indices from which we estimate the velocity
left_indices  = NaN(1,length(middle_indices));
right_indices = NaN(1,length(middle_indices));

%NOTE: This might be faster if we only tested unmatched values ...
%NOTE: We could also switch the matched logic to avoid the negation ...
for iShift = half_scale:scale_minus_1
    
    %We grab indices that are the appropriate distance from the current
    %value. If we have not yet found a bound on the given side, and the
    %index is valid, we keep it.
    left_indices_temp  = middle_indices - iShift;
    right_indices_temp = middle_indices + iShift;
    
    is_good_left       = is_good_value_mask(left_indices_temp);
    is_good_right      = is_good_value_mask(right_indices_temp);
    
    %Use if the value is good and not set ..
    use_left  = ~matched_left_mask  & is_good_left;
    use_right = ~matched_right_mask & is_good_right;
    
    left_indices(use_left)   = left_indices_temp(use_left);
    right_indices(use_right) = right_indices_temp(use_right);
    
    matched_left_mask(use_left)   = true;
    matched_right_mask(use_right) = true;
end

left_indices   = left_indices    - half_scale; %Remove the offset ...
right_indices  = right_indices   - half_scale;
middle_indices = middle_indices  - half_scale;

%Filter down to usable values
%NOTE: Some values may not have valid points on one or both sides ...
keep_mask = ~isnan(left_indices) & ~isnan(right_indices);

left_indices   = left_indices(keep_mask);
right_indices  = right_indices(keep_mask);
middle_indices = middle_indices(keep_mask);

%--------------------------------------------------------------------------

dX = x_mean(right_indices) - x_mean(left_indices);
dY = y_mean(right_indices) - y_mean(left_indices);
distance = sqrt(dX.^ 2 + dY.^ 2);
time     = (right_indices - left_indices)./ fps;
speed(middle_indices) = distance./time;

direction(middle_indices)   = pointAngle(right_indices) - pointAngle(left_indices);
direction(direction < -180) = direction(direction < -180) + 360;
direction(direction > 180)  = direction(direction > 180)  - 360;
direction = direction./fps;

% Sign the speed.
%--------------------------------------------------------
motionDirection = atan2(dY, dX) * (180 / pi);
bodyDirection(middle_indices) = motionDirection - bodyAngle(left_indices);

bodyDirection(bodyDirection < -180) = bodyDirection(bodyDirection < -180) + 360;
bodyDirection(bodyDirection > 180)  = bodyDirection(bodyDirection > 180)  - 360;

speed(abs(bodyDirection) > 90) = -speed(abs(bodyDirection) > 90);

% Sign the direction for dorsal/ventral locomtoion.
if ventralMode < 2 % + = dorsal direction
    direction = -direction;
end

% Organize the velocity.
velocity.speed     = speed;
velocity.direction = direction;
end




















%{
This is the original code, although I added some comments ...
function velocity = h__computeVelocity2(x, y, bodyAngle, pointsI, fps, scale, ventralMode)

%What is scale ????

% The scale must be odd.
scale = scale * fps;

%???? 
if rem(floor(scale), 2)
    scale = floor(scale);
elseif rem(ceil(scale), 2)
    scale = ceil(scale);
else
    scale = round(scale + 1);
end

% Do we have enough coordinates?
speed     = nan(1, size(x, 2));
direction = nan(1, length(speed));
body_direction = nan(1,length(speed));
if scale > size(x, 2)
    velocity.speed     = speed;
    velocity.direction = direction;
    return;
end

% Compute the body part direction.
diffX      = nanmean(diff(x(pointsI,:), 1, 1), 1);
diffY      = nanmean(diff(y(pointsI,:), 1, 1), 1);
pointAngle = atan2(diffY, diffX) * (180 / pi);

% Compute the coordinates.
x = mean(x(pointsI,:), 1);
y = mean(y(pointsI,:), 1);

% Compute the speed using back/front nearest neighbors bounded at twice the scale.
scaleMinus1 = scale - 1;
halfScale   = scaleMinus1 / 2;
diff1       = 1;
diff2       = scale;

for i = (1 + halfScale):(length(speed) - halfScale)
    
    % Advance the indices for the derivative.
    newDiff1 = i - halfScale;
    if ~isnan(x(newDiff1))
        diff1 = newDiff1;
    elseif i - diff1 >= scale
        diff1 = i - scale + 1;
    end
    newDiff2 = i + halfScale;
    if ~isnan(x(newDiff2)) || newDiff2 > diff2
        diff2 = newDiff2;
    end
    
    % Find usable indices for the derivative.
    while isnan(x(diff1)) && diff1 > 1 && i - diff1 < scaleMinus1
        diff1 = diff1 - 1;
    end
    while isnan(x(diff2)) && diff2 < length(speed) && diff2 - i < scaleMinus1
        diff2 = diff2 + 1;
    end
    
    % Compute the speed.
    if ~isnan(x(diff1)) && ~isnan(x(diff2))
        diffX    = x(diff2) - x(diff1);
        diffY    = y(diff2) - y(diff1);
        distance = sqrt(diffX ^ 2 + diffY ^ 2);
        time     = (diff2 - diff1) / fps;
        speed(i) = distance / time;
        
        % Compute the direction.
        direction(i) = pointAngle(diff2) - pointAngle(diff1);
        if direction(i) < -180
            direction(i) = direction(i) + 360;
        elseif direction(i) > 180
            direction(i) = direction(i) - 360;
        end
        direction(i) = direction(i) / fps;
        
        % Sign the speed.
        %--------------------------------------------------------
        motionDirection = atan2(diffY, diffX) * (180 / pi);
        bodyDirection   = motionDirection - bodyAngle(diff1);
        
        if bodyDirection < -180
            bodyDirection = bodyDirection + 360;
        elseif bodyDirection > 180
            bodyDirection = bodyDirection - 360;
        end
        if abs(bodyDirection) > 90
            speed(i) = -speed(i);
        end
        body_direction(i) = bodyDirection;
    end
end

% Sign the direction for dorsal/ventral locomtoion.
if ventralMode < 2 % + = dorsal direction
    direction = -direction;
end

% Organize the velocity.
velocity.speed     = speed;
velocity.direction = direction;
velocity.body_direction = body_direction;
end
%}