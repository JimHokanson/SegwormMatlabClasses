function wormPathCurvature(obj, x, y, fps, ventral_mode)
%
%   Compute the worm path curvature (angle/distance).
%
%   seg_worm.feature_helpers.path.wormPathCurvature
%
%   Old Name: wormPathCurvature.m
%
%   Inputs:
%   =======================================================================
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
%   =======================================================================
%       curvature - the worm path curvature (the angle between every 3 
%                   subsequent locations at the given scale, divided by the
%                   distance traveled between these 3 subsequent locations)
%
%   Nature Methods Description
%   =======================================================================
%   Curvature. 
%   -----------------------------
%   The path curvature is defined as the angle, in radians, of the worm’s
%   path divided by the distance it traveled in microns. The curvature is
%   signed to provide the path’s dorsal-ventral orientation. When the
%   worm’s path curves in the direction of its ventral side, the curvature
%   is signed negatively.
% 
%   The worm’s location is defined as the centroid of its body, with the
%   head and tail removed (points 9-41). We remove the head and tail
%   because their movement can cause large displacements in the worm’s
%   centroid. 
%
%   For each frame wherein the worm’s location is known, we search for a
%   start frame 1/4 of a second before and an end frame 1/4 second after to
%   delineate the worm’s instantaneous path. If the worm’s location is not
%   known within either the start or end frame, we extend the search for a
%   known location up to 1/2 second in either direction. If the worm’s
%   location is still missing at either the start or end, the path
%   curvature is marked unknown at this point.
% 
%   With three usable frames, we have an approximation of the start,
%   middle, and end for the worm’s instantaneous path curvature. We use the
%   difference in tangent angles between the middle to the end and between
%   the start to the middle. The distance is measured as the integral of
%   the distance traveled, per frame, between the start and end frames.
%   When a frame is missing, the distance is interpolated using the next
%   available segmented frame. The instantaneous path curvature is then
%   computed as the angle divided by the distance. This path curvature is
%   signed negatively if the angle curves in the direction of the worm’s
%   ventral side.
%
%   See Also:
%   

BODY_DIFF = 0.5; %Minimum window over which to calculate the velocity

% Initialize the body parts.
bodyI = 45:-1:5; %This does not match the description ...

% Compute the tail-to-head direction.
diffX     = nanmean(diff(x(bodyI,:), 1, 1), 1);
diffY     = nanmean(diff(y(bodyI,:), 1, 1), 1);
avg_body_angles_d = atan2(diffY, diffX) * (180 / pi);

velocity = seg_worm.features.helpers.computeVelocity(...
    x, y, avg_body_angles_d, bodyI, fps, BODY_DIFF, ventral_mode);

% Compute the path curvature.
obj.curvature = h__computeCurvature(velocity.speed, velocity.motion_direction, BODY_DIFF, fps);

end


%% Compute the worm path curvature.
function curvature = h__computeCurvature(speed, motion_direction, window_width, fps)
%
%
%   Inputs
%   ==============================================================
%   speed : [1 x n_frames]
%   motion_direction : [1 x n_frames] (Units: degrees)  TODO: Add description
%   window_width : (scalar) (Units: s)
%   fps : (scalar)
%
%   Outputs
%   ================================================
%   curvature

% The frame scale must be odd.
frame_scale = seg_worm.features.helpers.getWindowWidthAsInteger(window_width,fps);
half_frame_scale = (frame_scale - 1) / 2;

% Compute the angle differentials and distances.
speed = abs(speed);

diff_motion = NaN(size(speed));
right_max_I = length(diff_motion) - frame_scale + 1;
diff_motion(1:right_max_I) = motion_direction(frame_scale:end) - motion_direction(1:right_max_I);

diff_motion(diff_motion >= 180) = diff_motion(diff_motion >= 180) - 360;
diff_motion(diff_motion <= -180) = diff_motion(diff_motion <= -180) + 360;

distanceI           = (half_frame_scale + 1):(length(speed) - frame_scale);
distance            = NaN(size(speed));
distance(distanceI) = ((speed(distanceI) + speed(distanceI + frame_scale)) * window_width) / 2;

% Wrap the direction.


% Compute the worm path curvature.
distance(distance < 1) = NaN;
curvature = (diff_motion ./ distance) * (pi / 180);
end



