function getDirections
%
%
%   seg_worm.feature_helpers.posture.getDirections

% Head and tail centroids
headCentroid = mean(skCoords(1:round(1/6*NUMBER_OF_POINTS),:));
tailCentroid = mean(skCoords(round(5/6*NUMBER_OF_POINTS)+1:NUMBER_OF_POINTS,:));

% Compute tail direction
tailToHeadDirectionFrame   = atan2(headCentroid(1, 2) - tailCentroid(1,2), headCentroid(1,1) - tailCentroid(1,1));
tailToHeadDirection(frame) = tailToHeadDirectionFrame * 180/pi;

% Compute head and tail direction
headEnd           = 1:round(1/18*NUMBER_OF_POINTS);
headBegin         = round(1/6*NUMBER_OF_POINTS)+1 - headEnd;
headBegin         = fliplr(headBegin);
headEndCentroid   = mean(skCoords(headEnd,:));
headBeginCentroid = mean(skCoords(headBegin,:));

% Tail
tailEnd           = round(17/18*NUMBER_OF_POINTS) + 1:NUMBER_OF_POINTS;
tailBegin         = round(5/6*NUMBER_OF_POINTS) + 1:round(16/18*NUMBER_OF_POINTS);
tailEndCentroid   = mean(skCoords(tailEnd,:));
tailBeginCentroid = mean(skCoords(tailBegin,:));

% Direction for head
headDirectionFrame   = atan2(headEndCentroid(2) - headBeginCentroid(2), headEndCentroid(1) - headBeginCentroid(1));
headDirection(frame) = headDirectionFrame * 180/pi;

% Direction for tail
tailDirectionFrame   = atan2(tailEndCentroid(2) - tailBeginCentroid(2), tailEndCentroid(1) - tailBeginCentroid(1));
tailDirection(frame) = tailDirectionFrame * 180/pi;





% featureData.tailToHeadDirection = tailToHeadDirection;
% featureData.headDirection       = headDirection;
% featureData.tailDirection       = tailDirection;

% tailToHeadDirection = featureData.tailToHeadDirection;
% headPosDirection    = featureData.headDirection;
% tailPosDirection    = featureData.tailDirection;


% postureDirections = struct( ...
%     'tail2head',    tailToHeadDirection, ...
%     'head',         headPosDirection, ...
%     'tail',         tailPosDirection);


worm.posture.directions.tail2head
worm.posture.directions.head
worm.posture.directions.tail