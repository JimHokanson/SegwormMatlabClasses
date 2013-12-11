function locomotion = getLocomotionFeatures(nw)
%
%   locomotion = seg_worm.feature_calculator.getLocomotionFeatures(nw)
%
%   Inputs
%   =======================================================================
%   nw : seg_worm.normalized_worm
%
%
%   See Also:
%   seg_worm.feature_helpers.computeVelocity
%
%
%
%
%
%   TODO - the velocities are not the same ...

%ventralMode:
%0 - unknown
%1 - clockwise
%2 - anticlockwise

FPS = 25.8398;
VENTRAL_MODE = 0;  %??? How is this computed ????
%NOTE: midbody speed was opposite, which meant ventral mode was 2
%TODO: still need to get this from somewhere ...
%
%
%   Used in:
%   seg_worm.feature_helpers.locomotion.getWormVelocity
%

%Velocity - DONE
%--------------------------------------------------------------------------
locomotion.velocity = seg_worm.feature_helpers.locomotion.getWormVelocity(nw.x,nw.y,FPS,VENTRAL_MODE);


%Motion - DONE
%--------------------------------------------------------------------------
midbody_speed     = locomotion.velocity.midbody.speed;
locomotion.motion = seg_worm.feature_helpers.locomotion.getWormMotionCodes(midbody_speed, nw.lengths, FPS);


%Bends - DONE, needs documentation though ...
%--------------------------------------------------------------------------
motion_mode = locomotion.motion.mode;
is_paused   = motion_mode == 0;

%This is slow because of the significant # of ffts being computed ...
locomotion.bends = seg_worm.feature_helpers.locomotion.getLocomotionBends(...
    nw.angles, is_paused, nw.is_segmented, FPS);


locomotion.bends.foraging = ...
    seg_worm.feature_helpers.locomotion.getForaging(...
    nw.is_segmented,nw.x,nw.y,VENTRAL_MODE);


%Turns - still working on this ...
%--------------------------------------------------------------------------
midbody_distance = abs(midbody_speed/FPS);
is_stage_movement = nw.segmentation_status == 'm';
[omegas,upsilons] = ...
    seg_worm.feature_helpers.locomotion.getOmegaAndUpsilonTurns(...
    nw.angles,is_stage_movement,midbody_distance,nw.x,nw.y,FPS);

locomotion.turns.omegas   = omegas;
locomotion.turns.upsilons = upsilons;


keyboard










end