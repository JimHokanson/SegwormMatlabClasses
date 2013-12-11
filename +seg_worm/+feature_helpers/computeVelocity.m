function velocity = computeVelocity(sx, sy, avg_body_angle_d, pointsI, fps, scale_time, ventral_mode)
%computeVelocity  
%
%   seg_worm.feature_helpers.computeVelocity(x, y, bodyAngle, pointsI, fps, scale_time, ventralMode)
%
%   This function is needed by: 
%       seg_worm.feature_helpers.locomotion.getWormVelocities
%       
%
%   Old Name: 
%   - part of wormVelocity.m
%
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
%   Velocity. 
%   ---------------------------------------------
%   The worm’s velocity is measured at the tip of the head and
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

n_frames = size(sx,2);

%We need to go from a time over which to compute the velocity to a # of
%samples. The # of samples should be odd.
scale_samples_final = h__computeCalculationWidth(scale_time,fps);

% Do we have enough coordinates?
if scale_samples_final > n_frames
    velocity.speed     = nan(1,n_frames);
    velocity.direction = nan(1,n_frames);
    return;
end


%------------------------------------------------------------------
%Compute the indices that we will use for computing the velocity. We
%calculate the velocity roughly centered on each sample, but with a
%considerable width that smooths the velocity.
good_frames_mask = ~isnan(avg_body_angle_d);
[keep_mask,left_I,right_I] = h__getVelocityIndices(n_frames,scale_samples_final,good_frames_mask);

%Compute speed
%------------------------------------------------------------------
x_mean = mean(sx(pointsI,:), 1);
y_mean = mean(sy(pointsI,:), 1);

dX  = x_mean(right_I) - x_mean(left_I);
dY  = y_mean(right_I) - y_mean(left_I);

distance = sqrt(dX.^2 + dY.^2);
time     = (right_I - left_I)./ fps;

speed    = NaN(1,n_frames);
speed(keep_mask) = distance./time;

%Compute angular speed
%--------------------------------------------------------
angular_speed = NaN(1,n_frames); %Formally known as direction :/
angular_speed(keep_mask) = h__computeAngularSpeed(sx,sy,pointsI,left_I,right_I,ventral_mode,fps);

% Sign the speed.
%--------------------------------------------------------
%We want to know how the worm's movement direction compares to the average
%angle it had (apparently at the start)
motion_direction = atan2(dY, dX) * (180 / pi);

%This recenters the definition, as we are considered really with the
%change, not with the actual value
body_direction = NaN(1,n_frames);
body_direction(keep_mask) = motion_direction - avg_body_angle_d(left_I);

body_direction(body_direction < -180) = body_direction(body_direction < -180) + 360;
body_direction(body_direction > 180)  = body_direction(body_direction > 180)  - 360;

speed(abs(body_direction) > 90) = -speed(abs(body_direction) > 90);

% Organize the velocity.
%-----------------------------------------------------------
velocity.speed     = speed;
velocity.direction = angular_speed;

end

function scale_samples_final = h__computeCalculationWidth(scale_time,fps)

%The scale is used to determine over how many samples the velocity is
%calculated.

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

end

function angular_speed = h__computeAngularSpeed(sx,sy,pointsI,left_I,right_I,ventral_mode,fps)
%
%
%   
%


% Compute the body part direction.
avg_dX = nanmean(diff(sx(pointsI,:), 1, 1), 1);
avg_dY = nanmean(diff(sy(pointsI,:), 1, 1), 1);

point_angle_d = atan2(avg_dY, avg_dX) * (180 / pi);

angular_speed = point_angle_d(right_I) - point_angle_d(left_I);

%Correct any jumps that result during the subtraction process
%i.e. 1 - 359 ~= -358
angular_speed(angular_speed < -180) = angular_speed(angular_speed < -180) + 360;
angular_speed(angular_speed > 180)  = angular_speed(angular_speed > 180)  - 360;

angular_speed = angular_speed./fps;

% Sign the direction for dorsal/ventral locomtoion.
if ventral_mode < 2 % + = dorsal direction
    angular_speed = -angular_speed;
end


end

function [keep_mask,left_I,right_I] = h__getVelocityIndices(n_frames,scale_samples_final,good_frames_mask)
%
%
%
% Compute the speed using back/front nearest neighbors bounded at twice the scale.
%
%
%   Outputs
%   =======================================================================
%   keep_mask : [1 x n_frames], this is used to indicate which original
%               frames have valid velocity values, and which don't. 
%               NOTE: sum(keep_mask) == n_valid_velocity_values
%   left_I    : [1 x n_valid_velocity_values], for a given sample, this
%               indicates the index to the left of (less than) the sample
%               that should be used to calculate the velocity
%   right_I   : [1 x n_valid_velocity_values]


scale_minus_1 = scale_samples_final - 1;
half_scale    = scale_minus_1 / 2; %NOTE: Since the scale is odd, the half
%scale will be even, because we subtract 1

start_index = 1+half_scale; %First frame for which we can assign a valid velocity
end_index   = n_frames-half_scale;

%These are the indices we will use to compute the velocity. We add
%a half scale here to avoid boundary issues. We'll subtract it out later.
%See below for more details
middle_I     = (start_index:end_index) + half_scale;
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

unmatched_left_mask  = true(1,length(middle_I));
unmatched_right_mask = true(1,length(middle_I));

%This tells us whether each value is useable or not for velocity
%Better to do this out of the loop.
%For real indices (frames 1:n_frames), we rely on whether or not the
%mean position is NaN, for fake padding frames they can never be good so we
%set them to be false
is_good_value_mask = [false(1,half_scale) good_frames_mask false(1,half_scale)];

%These will be the final indices from which we estimate the velocity.
%i.e. delta_position(I) = (position(right_indices(I)) - position(left_indices(I))
left_I  = NaN(1,length(middle_I));
right_I = NaN(1,length(middle_I));

%Instead of looping over each centered velocity, we loop over each possible
%shift. A shift is valid if the frame of the shift is good, and we have yet
%to encounter a match for that centered index
for iShift = half_scale:scale_minus_1
    
    %We grab indices that are the appropriate distance from the current
    %value. If we have not yet found a bound on the given side, and the
    %index is valid, we keep it.
    left_indices_temp  = middle_I - iShift;
    right_indices_temp = middle_I + iShift;
    
    is_good_left_mask  = is_good_value_mask(left_indices_temp);
    is_good_right_mask = is_good_value_mask(right_indices_temp);
    
    use_left_mask      = unmatched_left_mask  & is_good_left_mask;
    use_right_mask     = unmatched_right_mask & is_good_right_mask;
    
    left_I(use_left_mask)   = left_indices_temp(use_left_mask);
    right_I(use_right_mask) = right_indices_temp(use_right_mask);
    
    unmatched_left_mask(use_left_mask)   = false;
    unmatched_right_mask(use_right_mask) = false;
end

left_I   = left_I    - half_scale; %Remove the offset ...
right_I  = right_I   - half_scale;
middle_I = middle_I  - half_scale;

%Filter down to usable values, in which both left and right are defined
valid_indices_mask = ~isnan(left_I) & ~isnan(right_I);
left_I    = left_I(valid_indices_mask);
right_I   = right_I(valid_indices_mask);
middle_I  = middle_I(valid_indices_mask);

keep_mask = false(1,n_frames);
keep_mask(middle_I) = true;

end
