%wormTouchFrames 
%   - drop and stage code switched ...
%   - last frame error is off by 1 (the bit at the end)
%     because 'i' doesn't advance like it does in the loop
%
%bends - indexing was incorrect for posture
%
%findEvent
%   sum data thresholding not implemented correctly
%   removeEvents using || should be |
%   off by 1 on start and end when event encompasses start or end