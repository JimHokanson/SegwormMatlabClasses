function locomotion = getLocomotionFeatures(nw)
%
%   locomotion = seg_worm.feature_calculator.getLocomotionFeatures(nw)
%
%   UNFINISHED STATUS
%   --------------------------
%   - bends not started
%   - turns started, but I need to incorporate new event code

%{
locomotion.bends.foraging.amplitude
locomotion.bends.foraging.angleSpeed

locomotion.bends.head.amplitude
locomotion.bends.midbody.amplitude
locomotion.bends.tail.amplitude
locomotion.bends.head.frequency
locomotion.bends.midbody.frequency
locomotion.bends.tail.frequency

locomotion.turns.omegas
locomotion.turns.upsilons

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


%Motion
%--------------------------------------------------------------------------
midbody_speed = locomotion.velocity.midbody.speed;
lengths       = nw.lengths;
%Almost done, needs a little more refactoring ...
locomotion.motion   = seg_worm.feature_helpers.locomotion.getWormMotionCodes(midbody_speed, lengths, FPS);

keyboard

%Bends
%--------------------------------------------------------------------------


%Turns
%--------------------------------------------------------------------------
[omegas,upsilons] = seg_worm.feature_helpers.locomotion.getOmegaAndUpsilonTurns(nw,1);

locomotion.turns.omegas   = omegas;
locomotion.turns.upsilons = upsilons;












end