function velocity = getWormVelocity(x, y, fps, ventralMode)
%getWormVelocity   Compute the worm velocity (speed & direction) at the
%head-tip/head/midbody/tail/tail-tip
%
%   velocity = seg_worm.feature_helpers.locomotion.getWormVelocity(x, y, fps, ventralMode)
%
%   Old Name: wormVelocity.m
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

SI        = seg_worm.skeleton_indices;
TIP_DIFF  = 0.25;
BODY_DIFF = 0.5;

% Compute the tail-to-head direction.
bodyI     = SI.BODY_INDICES;

diffX     = nanmean(diff(x(bodyI,:), 1, 1), 1);
diffY     = nanmean(diff(y(bodyI,:), 1, 1), 1);
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
   velocity.(cur_field_name)  = seg_worm.feature_helpers.computeVelocity(x, y, avg_body_angles_d, cur_indices, fps, cur_scale, ventralMode);
end

end



