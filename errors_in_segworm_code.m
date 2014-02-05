%wormTouchFrames 
%--------------------------------------------------------------------
%
%   - drop and stage code switched ...
%   - last frame error is off by 1 (the bit at the end)
%     because 'i' doesn't advance like it does in the loop
%
%bends - TODO: Clarify which function is being used
%--------------------------------------------------------------------
%- indexing was incorrect for posture
%
%
%findEvent
%--------------------------------------------------------------------
%   - sum data thresholding not implemented correctly
%   - event indices are 0 based, not 1 based 
%
%getAmpWavelength
%--------------------------------------------------------------------   
%   - power instead of magnitude is used for comparison
%   - primary and secondary wavelength may be switched ...
%   - primary and secondary both capped? - drop secondary in that case?
%
%
%seg_worm.feature_helpers.computeVelocity 
%--------------------------------------------------------------------
%- description in supplemental doesn't match reality ...
%
%
%seg_worm.feature_helpers.locomotion.getForaging
%--------------------------------------------------------------------
%
%   is the speed calculated correctly? Multiplying by fps???
%   I'm pretty sure it isn't correct
%
%
%seg_worm.feature_helpers.locomotion.getOmegaAndUpsilonTurns
%--------------------------------------------------------------------
%
%   Mismatch between description and cutoffs actually used for finding
%   possible event frames.
%
%
%seg_worm.feature_helpers.path.wormPathCurvature
%--------------------------------------------------------------------
% - indices used in body angle doesn't match description
% - NOTE: There is a comment about not using the ends because of noise,
%   but they are in seg_worm.feature_helpers.locomotion.getWormVelocity
%
%
%removePartialEvents.m
%--------------------------------------------------------------------------
%- indexing for the end event is off by 1
%
%worm2StatsInfo.m
%--------------------------------------------------------------------------
%- description of z-score doesn't match reality