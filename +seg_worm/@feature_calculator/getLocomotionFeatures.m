function locomotion = getLocomotionFeatures(nw)
%
%   locomotion = seg_worm.feature_calculator.getLocomotionFeatures(nw)
%
%   Unfinished  Unfinished  Unfinished  Unfinished

%{
locomotion.bends.foraging.amplitude
locomotion.bends.foraging.angleSpeed

locomotion.bends.head.amplitude
locomotion.bends.midbody.amplitude
locomotion.bends.tail.amplitude
locomotion.bends.head.frequency
locomotion.bends.midbody.frequency
locomotion.bends.tail.frequency

%}


%ventralMode:
%0 - unknown
%1 - clockwise
%2 - anticlockwise

FPS = 20;
VENTRAL_MODE = 1;  %??? How is this computed ????

%Velocity - DONE
%--------------------------------------------------------------------------
locomotion.velocity = seg_worm.feature_helpers.locomotion.getWormVelocity(nw.x,nw.y,FPS,VENTRAL_MODE);


%Motion - DONE
%--------------------------------------------------------------------------
midbody_speed     = locomotion.velocity.midbody.speed;
locomotion.motion = seg_worm.feature_helpers.locomotion.getWormMotionCodes(midbody_speed, nw.lengths, FPS);


%Bends - done, needs documentation though ...
%--------------------------------------------------------------------------
motion_mode = locomotion.motion.mode;
is_paused   = motion_mode == 0;
locomotion.bends = seg_worm.feature_helpers.locomotion.getLocomotionBends(...
    nw.angles, is_paused, nw.is_segmented);

locomotion.bends.foraging = ...
    seg_worm.feature_helpers.locomotion.getForaging(...
    nw.is_segmented,nw.x,nw.y,VENTRAL_MODE);

keyboard

%Turns - still working on this ...
%--------------------------------------------------------------------------
[omegas,upsilons] = seg_worm.feature_helpers.locomotion.getOmegaAndUpsilonTurns(nw,1);

locomotion.turns.omegas   = omegas;
locomotion.turns.upsilons = upsilons;












end