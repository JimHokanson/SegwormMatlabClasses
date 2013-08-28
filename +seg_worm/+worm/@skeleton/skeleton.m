classdef skeleton < handle
    %
    %   Class:
    %   seg_worm.worm.skeleton
    
    properties
        
        %SKELETON
        %------------------------------------------------------------------
        pixels %(sPixels)  - the worm's continuous skeleton oriented from head to tail
        touch_I %(sTouchI) - the paired pairs of indices marking, clockwise, the
        % start and end of the touching skeleton points
        % Note: if the worm isn't coiled, this value is empty.
        inner_I %(sInI) - the paired indices marking, clockwise, the start and
        % end of the inner skeleton points
        % Note: if the worm isn't coiled, this value is empty.
        outer_I %(sOutI) - the paired indices marking, clockwise, the start and
        % end of the outer skeleton points
        % Note: if the worm isn't coiled, this value is empty.
        inner_outer_I %(sInOutI) - the pairs of indices marking, from head to tail,
        % the start and end of the dual inner/outer skeleton points
        % Note: if the worm isn't coiled, this value is empty.
        angles %(sAngles) - the skeleton's angles (curvature) per point
        % Note 1: NaNs indicate the end pieces where there is
        % insufficient information to compute the angle
        % Note 2: positive skeleton angles bulge towards the side
        % clockwise from the worm's head (unless the worm is flipped)
        lengths %(sLength) - the skeleton's (worm's) chain-code pixel length
        chain_code_lengths %(sCCLengths) - the skeleton's (worm's) chain-code pixel length, from
        % head to tail, up to each skeleton point
        % Note: this is a more accurate representation of
        % locations along the worm's skeleton than pixel indices
        widths %(sWidths) - the contour's (worm's) widths, from head to tail, at
        % each skeleton point
    end
    
    methods
    end
    
end

