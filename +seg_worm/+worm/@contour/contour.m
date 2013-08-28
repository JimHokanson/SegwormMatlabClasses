classdef contour < handle
    %
    %   Class:
    %   seg_worm.worm.contour
    
    properties
        %CONTOUR
        %------------------------------------------------------------------
        pixels %(cPixels) the worm's circularly continuous contour pixels,
        %ordered clockwise
        touch_points_I  %(cTouchI) the paired pairs of indices marking, clockwise, the
        % start and end of the touching contour points Note: if the worm isn't
        % coiled, this value is empty.
        inner_I %(cInI) the paired indices marking, clockwise, the start and
        % end of the inner contour points
        outer_I %(cOutI) the paired indices marking, clockwise, the start and
        % end of the outer contour points
        angles  %(cAngles) the contour's angles (curvature) at each index
        head_I  %(cHeadI)  the contour index for the worm's head
        tail_I  %(cTailI)  the contour index for the worm's tail
        chain_code_lengths %(cCCLengths) - the contour's circular chain-code pixel length, from
        % its vector's start to end, up to each contour point Note: this is a more
        % accurate representation of locations along the worm's contour than pixel
        % indices
    end
    
    methods
    end
    
end

