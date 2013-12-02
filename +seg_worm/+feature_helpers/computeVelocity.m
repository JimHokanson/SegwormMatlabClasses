function velocity = computeVelocity(x, y, avg_body_angle_d, pointsI, fps, scale_time, ventralMode)
%computeVelocity  
%
%   seg_worm.feature_helpers.computeVelocity(x, y, bodyAngle, pointsI, fps, scale_time, ventralMode)
%
%   This function is needed by: 
%       seg_worm.feature_helpers.locomotion.getWormVelocities
%
%   Original MRC File: wormVelocity.m
%
%
%   Improvements to Make:
%   =======================================================================
%   1: 
%
%   Inputs
%   =======================================================================
%   x           : [49 x n_frames], skeleton x points
%   y           : [49 x n_frames], skeleton y points
%   avg_body_angle_d : [1  x n_frames] - average body angle in degrees,
%                   each angle is computed using neighboring pairs of
%                   points (each pair forms a vector), with the angle being 
%                   calculated relative to the reference frame of the system, 
%                   not another vector, the vector tip is closer to the
%                   woirm tail than the base of the vector
%   pointsI     : indices of points on the skeleton to use for the
%                 calculation
%   fps         : (scaler) frames per second
%   scale_time  : (scaler) (Units: s), see comment below on computing the velocity
%   ventralMode : (scaler) the ventral side mode:
%
%                     0 = the ventral side is unknown
%                     1 = the ventral side is clockwise
%                     2 = the ventral side is anticlockwise
%
%
%   Outputs
%   =======================================================================
%   velocity    : (struct)
%       .speed     - [1 x n_frames]
%       .direction - [1 x n_frames]
%
%   Computing the Velocity
%   =======================================================================
%   The velocity is computed not using the nearest values but values
%   that are separated by a sufficient time (scale_time). If the encountered 
%   values are not valid, the width of time is expanded up to a maximum of
%   2*scale_time (or technically, 1 scale time in either direction)
%
%   Thus, the scale_time is the time over which we compute the velocity. If
%   the scale_time is 1, the velocity at time t = 2 is computed using the
%   data from frames that are at times t = 1.5 and t = 2.5 . If the frames
%   at those times are bad we expand up to maximum bounds of t = 1 and t =
%   3. Note, it is fine to use lopsided values such as t = 1.3 and t = 2.5.
%  
%
%   Nature Methods Description
%   =======================================================================
%   Velocity. The worm’s velocity is measured at the tip of the head and
%   tail, at the head and tail themselves, and at the midbody. The velocity
%   is composed of two parts, speed and direction (expressed as an angular
%   speed) (Supplementary Fig. 4d). The velocity is signed negatively
%   whenever the respective body part moves towards the tail (as opposed to
%   the head).
%
%   The head and tail tips’ instantaneous velocity is measured at each
%   frame using a 1/4 second up to a 1/2 second window. For each frame, we
%   search for a start frame 1/4 of a second before and an end frame 1/4
%   second after to delineate the worm’s instantaneous path. If the worm’s
%   location is not known within either the start or end frame, we extend
%   the search for a known location up to 1/2 second in either direction.
%   If the worm’s location is still missing at either the start or end, the
%   velocity is marked unknown at this point. 
%
%   The speed is defined as the distance between the centroids of the start
%   and end frames (for the respective body parts) divided by the time
%   between both frames.
%
%   The direction is defined as the angle (between centroids) from the
%   start to the end frame, relative to the worm’s overall body angle (JAH:
%   ???, not sure if this is true), divided by the time between both
%   frames.
%
%   The worm’s overall body angle is defined as the mean orientation of the
%   angles, in the tail-to-head direction, between subsequent midbody
%   skeleton points. The body angle is used to sign the velocity. If the
%   head or tail tip’s start-to-end angle exceeds 90°, clockwise or
%   anticlockwise, relative to the overall worm body angle, the motion is
%   towards the tail. In this case both the speed and direction are
%   negatively signed (JAH: The direction is not signed based on the body
%   angle but rather on the dorsal-ventral definition).
%
%
%   The head, midbody, and tail velocity are computed
%   identically except they use a 1/2 second up to a 1 second window for
%   choosing their start and end frames.
%

n_frames = size(x,2);

speed         = nan(1,n_frames);
angular_speed = nan(1,n_frames); %Formally known as direction :/
bodyDirection = nan(1,n_frames);

% The scale must be odd.
%--------------------------------------------------------------------------
%To this end we take the value obtained and round it down and up, we choose
%the odd value. This is made tricky if it is an even integer to start, such
%that rounding up or down both result in an even integer. In this case we
%add 1
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
%--------------------------------------------------------------------------

% Do we have enough coordinates?
if scale_samples_final > n_frames
    velocity.speed     = speed;
    velocity.direction = angular_speed;
    return;
end

% Compute the body part direction.
dX = nanmean(diff(x(pointsI,:), 1, 1), 1);
dY = nanmean(diff(y(pointsI,:), 1, 1), 1);

point_angle_d = atan2(dY, dX) * (180 / pi);

% Compute the coordinates.
x_mean = mean(x(pointsI,:), 1);
y_mean = mean(y(pointsI,:), 1);

% Compute the speed using back/front nearest neighbors bounded at twice the scale.
%--------------------------------------------------------------------------
%
%   This is the tricky part ...
scale_minus_1 = scale_samples_final - 1;
half_scale    = scale_minus_1 / 2; %NOTE: Since the scale is odd, the half
%scale will be even, because we subtract 1

start_index = 1+half_scale; %First frame for which we can assign a valid velocity
end_index   = n_frames-half_scale;

%These are the indices we will use to compute the velocity. We add
%a half scale here to avoid boundary issues. We'll subtract it out later.
%See below for more details
middle_indices     = (start_index:end_index) + half_scale;
%   Our start_index frame can only have one valid start frame (frame 1)
%   afterwards it is invalid. In other words, if frame 1 is not good, we
%   can't check frame 0 or frame -1, or -2.
%
%   However, in general I'd prefer to avoid doing some check on the bounds
%   of the frames, i.e. for looking at starting frames, I'd like to avoid
%   checking if the frame value is 1 or greater.
%
%   To get around this we'll pad with bad values (which we'll never match)
%   then shift the final indices. In this way, we can check these "frames",
%   as they will have positive indices.
%
%   e.g.
%   scale = 5
%   half_scale = 2
%
%   This means the first frame in which we could possibly get a valid
%   velocity is frame 3, computed using frames 1 and 5
%
%   NaN NaN 1 2 3  <- true indices (frame numbers)
%   1   2   3 4 5  <- temporary indices
%
%   NOTE: Now if frame 1 is bad, we can go left by half_scale + 1 to temp
%   index 2 (frame 0) or even further to temp_index 1 (frame -1), we'll
%   never use those values however because below we state that the values
%   at those indices are bad (see is_good_value_mask)
%
%   

unmatched_left_mask  = true(1,length(middle_indices));
unmatched_right_mask = true(1,length(middle_indices));

%This tells us whether each value is useable or not for velocity
%Better to do this out of the loop.
%For real indices (frames 1:n_frames), we rely on whether or not the
%mean position is NaN, for fake padding frames they can never be good so we
%set them to be false
is_good_value_mask = [false(1,half_scale) ~isnan(x_mean) false(1,half_scale)];

%These will be the final indices from which we estimate the velocity.
%i.e. delta_position(I) = (position(right_indices(I)) - position(left_indices(I))
left_indices  = NaN(1,length(middle_indices));
right_indices = NaN(1,length(middle_indices));

%Instead of looping over each centered velocity, we loop over each possible
%shift. A shift is valid if the frame of the shift is good, and we have yet
%to encounter a match for that centered index
for iShift = half_scale:scale_minus_1
    
    %We grab indices that are the appropriate distance from the current
    %value. If we have not yet found a bound on the given side, and the
    %index is valid, we keep it.
    left_indices_temp  = middle_indices - iShift;
    right_indices_temp = middle_indices + iShift;
    
    is_good_left_mask  = is_good_value_mask(left_indices_temp);
    is_good_right_mask = is_good_value_mask(right_indices_temp);
    
    use_left_mask      = unmatched_left_mask  & is_good_left_mask;
    use_right_mask     = unmatched_right_mask & is_good_right_mask;
    
    left_indices(use_left_mask)   = left_indices_temp(use_left_mask);
    right_indices(use_right_mask) = right_indices_temp(use_right_mask);
    
    unmatched_left_mask(use_left_mask)   = false;
    unmatched_right_mask(use_right_mask) = false;
end

left_indices   = left_indices    - half_scale; %Remove the offset ...
right_indices  = right_indices   - half_scale;
middle_indices = middle_indices  - half_scale;

%Filter down to usable values, in which both left and right are defined
keep_mask      = ~isnan(left_indices) & ~isnan(right_indices);
left_indices   = left_indices(keep_mask);
right_indices  = right_indices(keep_mask);
middle_indices = middle_indices(keep_mask);

%Now onto the final calculations ...
%--------------------------------------------------------------------------
dX       = x_mean(right_indices) - x_mean(left_indices);
dY       = y_mean(right_indices) - y_mean(left_indices);
distance = sqrt(dX.^ 2 + dY.^ 2);
time     = (right_indices - left_indices)./ fps;

speed(middle_indices) = distance./time;

angular_speed(middle_indices)   = point_angle_d(right_indices) - point_angle_d(left_indices);
angular_speed(angular_speed < -180) = angular_speed(angular_speed < -180) + 360;
angular_speed(angular_speed > 180)  = angular_speed(angular_speed > 180)  - 360;

angular_speed = angular_speed./fps;

% Sign the direction for dorsal/ventral locomtoion.
if ventralMode < 2 % + = dorsal direction
    angular_speed = -angular_speed;
end


% Sign the speed.
%--------------------------------------------------------
%We want to know how the worm's movement direction compares to the average
%angle it had (apparently at the start)
motionDirection = atan2(dY, dX) * (180 / pi);

%This recenters the definition, as we are considered really with the
%change, not with the actual value
bodyDirection(middle_indices) = motionDirection - avg_body_angle_d(left_indices);

bodyDirection(bodyDirection < -180) = bodyDirection(bodyDirection < -180) + 360;
bodyDirection(bodyDirection > 180)  = bodyDirection(bodyDirection > 180)  - 360;

speed(abs(bodyDirection) > 90) = -speed(abs(bodyDirection) > 90);

% Organize the velocity.
%-----------------------------------------------------------
velocity.speed     = speed;
velocity.direction = angular_speed;
end










%==========================================================================
%                            THE OLD CODE
%==========================================================================





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

        %??????????? - why is this done ????
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