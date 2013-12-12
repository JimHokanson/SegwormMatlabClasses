function locomotion = getLocomotionFeatures(nw,FPS,VENTRAL_MODE)
%
%   seg_worm.feature_calculator.getLocomotionFeatures
%
%   Inputs
%   =======================================================================
%   nw : seg_worm.normalized_worm
%   FPS :
%   VENTRAL_MODE :
%
%   STATUS:
%   DONE
%
%   See Also:
%   seg_worm.feature_helpers.locomotion.getWormVelocity
%   seg_worm.feature_helpers.locomotion.getWormMotionCodes
%   seg_worm.feature_helpers.locomotion.getLocomotionBends
%   seg_worm.feature_helpers.locomotion.getForaging
%   seg_worm.feature_helpers.locomotion.getOmegaAndUpsilonTurns
%   seg_worm.feature_helpers.computeVelocity



%Velocity - DONE
%--------------------------------------------------------------------------
locomotion.velocity = seg_worm.feature_helpers.locomotion.getWormVelocity(nw.x,nw.y,FPS,VENTRAL_MODE);



%Motion - DONE
%--------------------------------------------------------------------------
midbody_speed     = locomotion.velocity.midbody.speed;
locomotion.motion = seg_worm.feature_helpers.locomotion.getWormMotionCodes(midbody_speed, nw.lengths, FPS);



%Crawling - part of bends - DONE
%--------------------------------------------------------------------------
motion_mode = locomotion.motion.mode;
is_paused   = motion_mode == 0;
%SLOW - many ffts being computed
locomotion.bends = seg_worm.feature_helpers.locomotion.getLocomotionBends(...
    nw.angles, is_paused, nw.is_segmented, FPS);



%Foraging - part of bends - DONE
%--------------------------------------------------------------------------
locomotion.bends.foraging = ...
    seg_worm.feature_helpers.locomotion.getForaging(...
    nw.x,nw.y,nw.is_segmented,VENTRAL_MODE,FPS);



%Turns - nearly done, could benefit from a bit more documentation
%--------------------------------------------------------------------------
midbody_distance  = abs(midbody_speed/FPS);
is_stage_movement = nw.segmentation_status == 'm';

[omegas,upsilons] = ...
    seg_worm.feature_helpers.locomotion.getOmegaAndUpsilonTurns(...
    nw.angles,is_stage_movement,midbody_distance,nw.x,nw.y,FPS);

locomotion.turns.omegas   = omegas;
locomotion.turns.upsilons = upsilons;

end