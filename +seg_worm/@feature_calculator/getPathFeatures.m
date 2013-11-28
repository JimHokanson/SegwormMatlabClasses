function path = getPathFeatures(nw)
%
%
%   seg_worm.feature_calculator.getPathFeatures(nw)
%
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


FPS = 20;

% % % %% The worm path.
% % % wormPath = struct( ...
% % %     'range',        pathRange, ...
% % %     'duration',     pathDuration, ...
% % %     'coordinates',  centroidCoordinates, ...
% % %     'curvature',    pathCurvature);


%Range
%---------------------------------------------------------------
path.range = seg_worm.feature_helpers.path.getRange(nw);




%Duration
%--------------------------------------------------------------------------
%       arena: [1x1 struct]
%        worm: [1x1 struct]
%        head: [1x1 struct]
%     midbody: [1x1 struct]
%        tail: [1x1 struct]

%{
worm.path.duration

%TODO: This information needs to come from the video I think ...
d.arena.height;
d.arena.width;
d.arena.min.x;
d.arena.min.y;
d.arena.max.x;
d.arena.max.y;

d.worm.indices
d.worm.times

d.head.indices
d.head.times

d.midbody.indices
d.midbody.times

d.tail.indices
d.tail.times

%}



path.duration = seg_worm.feature_helpers.path.getDurationInfo(nw.x, nw.y, nw.widths, FPS);


%Coordinates
%--------------------------------------------------------------------------
path.coordinates.x = mean(nw.contour_x);
path.coordinates.y = mean(nw.contour_y);




% Compute the normalized moments
pathX = outlineCentroid(:,1);
pathY = outlineCentroid(:,2);
%calculate centroid
pathCentroidX = nanmean(pathX);
pathCentroidY = nanmean(pathY);
distanceFromCenter = sqrt((pathX - pathCentroidX).^2 + (pathY - pathCentroidY).^2);
% Get the std and skewness of the distances
pathStd = nanstd(distanceFromCenter);
%pathSkewness = nanmean((distance/pathStd).^3);
featureData.pathStd = pathStd;

featureData.pathSkewness = skewness(distanceFromCenter, 0);
featureData.pathKurtosis = kurtosis(distanceFromCenter, 0);


centroidPathX = featureData.outlineCentroid(1,:);
centroidPathY = featureData.outlineCentroid(2,:);


%% Compute the path range.
pathRange = wormPathRange(centroidPathX, centroidPathY);

%% Compute the path duration.

% Compute the path scale.
headWidth = nanmean(headWidths);
midWidth  = nanmean(midbodyWidths);
tailWidth = nanmean(tailWidths);
meanWidth = (headWidth + midWidth + tailWidth) / 3;
pathScale = sqrt(2) / meanWidth;

% Compute the worm segments.
headI = 1;
tailI = wormSamples;
wormSegSize = round(tailI / 6);
headIs      = headI:(headI + wormSegSize - 1);
midbodyIs   = (headI + wormSegSize):(tailI - wormSegSize);
tailIs      = (tailI - wormSegSize + 1):tailI;

% Compute the skeleton points.
points = { ...
    headI:tailI, ...
    headIs, ...
    midbodyIs, ...
    tailIs};

% Compute the path duration and organize everything for the feature file.
[arena, durations] = wormPathTime(postureXSkeletons, postureYSkeletons, points, pathScale, fps);
pathDuration = struct( ...
    'arena',    arena, ...
    'worm',     durations(1), ...
    'head',     durations(2), ...
    'midbody',  durations(3), ...
    'tail',     durations(4));

