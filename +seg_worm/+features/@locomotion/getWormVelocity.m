function getWormVelocity(obj, sx, sy, fps, ventral_mode)
%getWormVelocity   Compute the worm velocity (speed & direction) at the
%head-tip/head/midbody/tail/tail-tip
%
%   seg_worm.feature_helpers.locomotion.getWormVelocity
%
%   Old Name: wormVelocity.m
%
%   Inputs:
%   =======================================================================
%       sx           : the worm skeleton's x-axis coordinates
%       sy           : the worm skeleton's y-axis coordinates
%       fps          : the frames/seconds
%       ventral_mode : the ventral side mode:
%
%                     0 = the ventral side is unknown
%                     1 = the ventral side is clockwise
%                     2 = the ventral side is anticlockwise
%
%   Outputs:
%   =======================================================================
%       velocity - the worm velocity; each field has subfields for the
%                  "speed" and "direction":
%
%                  headTip = the tip of the head (1/12 the worm at 0.25s)
%                  head    = the head (1/6 the worm at 0.5s)
%                  midbody = the midbody (2/6 the worm at 0.5s)
%                  tail    = the tail (1/6 the worm at 0.5s)
%                  tailTip = the tip of the tail (1/12 the worm at 0.25s)
%
%   Nature Methods Description
%   =======================================================================
%   Velocity. 
%   ----------------------------------
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
%   velocity is marked unknown at this point. The speed is defined as the
%   distance between the centroids of the start and end frames (for the
%   respective body parts) divided by the time between both frames. The
%   direction is defined as the angle (between centroids) from the start to
%   the end frame, relative to the worm’s overall body angle, divided by
%   the time between both frames. The worm’s overall body angle is defined
%   as the mean orientation of the angles, in the tail-to-head direction,
%   between subsequent midbody skeleton points. The body angle is used to
%   sign the velocity. If the head or tail tip’s start-to-end angle exceeds
%   90°, clockwise or anticlockwise, relative to the overall worm body
%   angle, the motion is towards the tail. In this case both the speed and
%   direction are negatively signed. The head, midbody, and tail velocity
%   are computed identically except they use a 1/2 second up to a 1 second
%   window for choosing their start and end frames.
%
%   See Also:
%   seg_worm.feature_helpers.computeVelocity

SI        = seg_worm.skeleton_indices;
TIP_DIFF  = 0.25;
BODY_DIFF = 0.5;

% Compute the tail-to-head direction.
%--------------------------------------------------------------------------
bodyI     = SI.BODY_INDICES(end:-1:1); %flip for angle calculations

%NOTE: This is different than in:
%seg_worm.feature_helpers.path.wormPathCurvature
diffX     = nanmean(diff(sx(bodyI,:), 1, 1), 1);
diffY     = nanmean(diff(sy(bodyI,:), 1, 1), 1);
avg_body_angles_d = atan2(diffY, diffX).*(180 / pi);

% Compute the velocity.
%--------------------------------------------------------------------------
FIELD_NAMES       = {'headTip'           'head'          'midbody'       'tail'          'tailTip'};
INDICES_LIST      = {SI.HEAD_TIP_INDICES SI.HEAD_INDICES SI.MID_INDICES SI.TAIL_INDICES SI.TAIL_TIP_INDICES};
TIME_SCALE_VALUES = [TIP_DIFF BODY_DIFF BODY_DIFF BODY_DIFF TIP_DIFF];

for iField = 1:length(FIELD_NAMES)
   cur_field_name = FIELD_NAMES{iField};
   cur_indices    = INDICES_LIST{iField};
   cur_scale      = TIME_SCALE_VALUES(iField);
   
   %NOTE: This was moved to a separate function because I found that this
   %same function was being used elsewhere in:
   %    seg_worm.feature_helpers.path.wormPathCurvature
   temp = seg_worm.feature_helpers.computeVelocity(...
       sx, sy, avg_body_angles_d, cur_indices, fps, cur_scale, ventral_mode);
   
   velocity.(cur_field_name) = struct('speed',temp.speed,'direction',temp.angular_speed);
end

obj.velocity = velocity;

end



