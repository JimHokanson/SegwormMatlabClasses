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


%Bends - still working on this ...
%--------------------------------------------------------------------------
% % % locomotion.bends.foraging.amplitude
% % % locomotion.bends.foraging.angleSpeed
% % % 
% % % locomotion.bends.head.amplitude
% % % locomotion.bends.midbody.amplitude
% % % locomotion.bends.tail.amplitude
% % % locomotion.bends.head.frequency
% % % locomotion.bends.midbody.frequency
% % % locomotion.bends.tail.frequency

%JAH: At this point
% - code needs to be simplified considerably ...
motion_mode = locomotion.motion.mode;
locomotion.bends = seg_worm.feature_helpers.locomotion.getLocomotionBends(...
    [], motion_mode, VENTRAL_MODE);

%This part is done ...
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