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
%   See Also:
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

%This is slow because of the significant # of ffts being computed.
locomotion.bends = seg_worm.feature_helpers.locomotion.getLocomotionBends(...
    nw.angles, is_paused, nw.is_segmented, FPS);



%Foraging - part of bends - needs documentation
%--------------------------------------------------------------------------
locomotion.bends.foraging = ...
    seg_worm.feature_helpers.locomotion.getForaging(...
    nw.is_segmented,nw.x,nw.y,VENTRAL_MODE);



%Turns - needs documentation
%--------------------------------------------------------------------------
midbody_distance = abs(midbody_speed/FPS);
is_stage_movement = nw.segmentation_status == 'm';
[omegas,upsilons] = ...
    seg_worm.feature_helpers.locomotion.getOmegaAndUpsilonTurns(...
    nw.angles,is_stage_movement,midbody_distance,nw.x,nw.y,FPS);

locomotion.turns.omegas   = omegas;
locomotion.turns.upsilons = upsilons;

end