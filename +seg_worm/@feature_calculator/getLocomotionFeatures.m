function locomotion = getLocomotionFeatures(nw)
%
%   locomotion = seg_worm.feature_calculator.getLocomotionFeatures(nw)
%

%{
locomotion.motion.forward
locomotion.motion.backward
locomotion.motion.paused
locomotion.motion.mode

locomotion.velocity.headTip.speed
locomotion.velocity.head.speed
locomotion.velocity.midbody.speed
locomotion.velocity.tail.speed
locomotion.velocity.tailTip.speed

locomotion.velocity.headTip.direction
locomotion.velocity.head.direction
locomotion.velocity.midbody.direction
locomotion.velocity.tail.direction
locomotion.velocity.tailTip.direction

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

%[velocity, motionEvents] = wormVelocity(postureXSkeletons, postureYSkeletons, fps, featureData.wormLength,  ventralMode);

%NOTE: I might break this up into multiple functions depending on what I
%see when I parse the function

%ventralMode:
%0 - unknown
%1 - clockwise
%2 - anticlockwise

FPS = 20;

ventralMode = 1;

locomotion.velocity = ...
    seg_worm.feature_helpers.locomotion.getWormVelocity(...
    nw.x,nw.y,FPS,ventralMode);

%,locomotion.motion]

end