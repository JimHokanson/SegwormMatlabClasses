function path = getPathFeatures(nw,FPS,VENTRAL_MODE)
%
%
%   seg_worm.feature_calculator.getPathFeatures(nw)
%
%worm.path.curvature - [1 x n]
%
%
% path. The path features. The path is represented by its “range”, “curvature”, and
% the dwelling “duration” for various body parts. Individual experiment files also
% contain the “x” and “y” “coordinates” of the contour’s centroid. Moreover, the
% individual experiment files present the “duration” as an “arena” with a “height”,
% “width”, and the “min” and “max” values for the “x” and “y” axes of the arena.
% The arena can be transformed to a matrix using the given height and width. The
% duration of the worm and body parts are represented as an array of “times” spent
% at the “indices” of the arena matrix
%
%
%   See also:
%   seg_worm.feature_helpers.path.getRange



%Range - DONE
%---------------------------------------------------------------
path.range = seg_worm.feature_helpers.path.getRange(nw.contour_x,nw.contour_y);


%Duration (aka Dwelling) - DONE
%--------------------------------------------------------------------------
path.duration = seg_worm.feature_helpers.path.getDurationInfo(nw.x, nw.y, nw.widths, FPS);


%Coordinates - DONE
%--------------------------------------------------------------------------
path.coordinates.x = mean(nw.contour_x);
path.coordinates.y = mean(nw.contour_y);


%Curvature - DONE
%--------------------------------------------------------------------------
path.curvature = seg_worm.feature_helpers.path.wormPathCurvature(nw.x,nw.y,FPS,VENTRAL_MODE);



