function path = getPathFeatures(nw)

%worm.path.range - [1 x n]
%worm.path.duration - struct
%
%         arena: [1x1 struct]
%                     height: 347
%                      width: 197
%                        min: [1x1 struct]
%                        max: [1x1 struct]        
%   
%        worm: [1x1 struct]
%                 indices: [2987x1 double]
%                   times: [2987x1 double]
%        head: [1x1 struct]
%               indices: [1924x1 double]
%               times: [1924x1 double]
%
%     midbody: [1x1 struct]
%           - same
%        tail: [1x1 struct]
%           - same
%
%   contours centroid coordinates
%worm.path.coordinates.x
%worm.path.coordinates.y
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