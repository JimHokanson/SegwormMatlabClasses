classdef movement_parser < handle
    %
    %   Class:
    %   seg_worm.stage.movement_parser
    %
    %   This will be the new home for frame parsing code ...
    %
    %   See Also:
    %   seg_worm.stage
    
    properties
       gOtsuThr
       gSmallThr   %median(gSmallDiffs) + 3 * std(gSmallDiffs);
       gSmallDiffs %frameDiffs(frameDiffs < gOtsuThr);
    end
    
    methods
    end
    
end

