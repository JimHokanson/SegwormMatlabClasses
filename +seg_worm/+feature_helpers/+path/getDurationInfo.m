function duration_struct = getDurationInfo(sx, sy, widths, fps)
%WORMPATHTIME Compute the time spent at each point along the worm path.
%
%   
%   duration_struct = seg_worm.feature_helpers.path.getDurationInfo(sx, sy, fps)
%
%   [ARENA TIMES] = WORMPATHTIME(SKELETON, POINTS, SCALE, FPS)
%
%   Inputs:
%       skeletonX - the worm skeleton's x coordinates per frame
%                   (x-coordinates x frames)
%       skeletonY - the worm skeleton's y coordinates per frame
%                   (y-coordinates x frames)
%       points    - the skeleton points (cell array) to use; each set of
%                   skeleton points delineates a separate path
%       scale     - the coordinate scale for the path
%                   Note: if the path time is integrated at the standard
%                   micron scale, the matrix will be too large and the
%                   worm's width will be too skinny. Therefore, I suggest
%                   using an a scale of 1/width to match the worm's width
%                   to a pixel edge; or, if you want to account for the
%                   long diagonal axis in the taxi-cab metric of pixels,
%                   use sqrt(2)/width to match the worm's width to a
%                   pixel's diagonal length
%       fps       - the frames/seconds
%
%   Output:
%       arena - a struct of the arena/path size with subfields:
%
%               height = the arena height
%                        (for the matrix of time spent at each point)
%               width  = the arena width
%                        (for the matrix of time spent at each point)
%
%               min:
%                  x = the path location of the arena's minimum x coordinate
%                  y = the path location of the arena's minimum y coordinate
%
%               max:
%                  x = the path location of the arena's maximum x coordinate
%                  y = the path location of the arena's maximum y coordinate
%
%       times - a struct(s) of the time(s) spent, per path, with subfields:
%
%               indices = the indices for the non-zero time points in the
%                         arena matrix
%               times   = the non-zero time point values (in seconds)
%                         corresponding to the arena matrix indices
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.


%Compute the scale
%--------------------------------------------------------------------------


% % % % headWidths    = featureData.widthsAtTips(1,:);
% % % % midbodyWidths = featureData.width;
% % % % tailWidths    = featureData.widthsAtTips(2,:);
% % % % 
% % % % headWidth = nanmean(headWidths);
% % % % midWidth  = nanmean(midbodyWidths);
% % % % tailWidth = nanmean(tailWidths);
% % % % meanWidth = (headWidth + midWidth + tailWidth) / 3;
% % % % scale = sqrt(2) / meanWidth;


%Compute the skeleton points
% % % 
% % % headI = 1;
% % % tailI = 49;
% % % wormSegSize = round(tailI / 6);
% % % headIs      = headI:(headI + wormSegSize - 1);
% % % midbodyIs   = (headI + wormSegSize):(tailI - wormSegSize);
% % % tailIs      = (tailI - wormSegSize + 1):tailI;

% % % %     headI:tailI, ...
% % % %     headIs, ...
% % % %     midbodyIs, ...
% % % %     tailIs};


SI = seg_worm.skeleton_indices;

% Compute the skeleton points.
s_points = {SI.ALL_INDICES SI.HEAD_INDICES SI.MID_INDICES SI.TAIL_INDICES};
n_points = length(s_points);

%??? - why scale this ????
mean_width = nanmean(mean(widths([s_points{2:4}],:),1),2);
scale = sqrt(2)/mean_width;
%NOTE: The old code omitted the widths at the neck and hips, I'm also doing
%that here, which is why we use the head, midbody, and tail indices,
%instead of just taking the mean of the widths. This distinction seems not
%functionally useful but I'll keep it for now ...


% The skeletons are empty.
if isempty(sx) || all(isnan(sx(:)))
    arena.height = NaN;
    arena.width = NaN;
    arena.min.x = NaN;
    arena.min.y = NaN;
    arena.max.x = NaN;
    arena.max.y = NaN;
    
    NAN_cell = repmat({NaN},1,n_points);
    
    durations = struct('indices',NAN_cell,'times',NAN_cell);
    duration_struct = h__buildOutput(arena,durations);
    return;
end

% Scale the skeleton.

sxs = round(sx*scale);
sys = round(sy*scale);

% Translate the skeleton to a zero origin.
xScaledMin = min(sxs(:));
xScaledMax = max(sxs(:));
yScaledMin = min(sys(:));
yScaledMax = max(sys(:));
sxs = sxs - xScaledMin + 1;
sys = sys - yScaledMin + 1;

% Compute the arena.
%--------------------------------------------------------------------------

% Construct the empty arena(s).
arenaSize = [yScaledMax - yScaledMin + 1, xScaledMax - xScaledMin + 1];

% Organize the arena size.
arena.height = arenaSize(1);
arena.width  = arenaSize(2);
arena.min.x  = min(sx(:));
arena.min.y  = min(sy(:));
arena.max.x  = max(sx(:));
arena.max.y  = max(sy(:));
%--------------------------------------------------------------------------


zeroArena = zeros(arenaSize);
arenas    = cell(length(s_points),1);
for iPoint = 1:n_points
    arenas{iPoint} = zeroArena;
end

n_frames = size(sxs,2);


%JAH TODO: I'm at this point rewriting the code ...
%
%NOTE: Instead of sorting, just do assignments
%
%   i.e.
%   
%   temp([1 3 5 5 3]) = 1;
%
%   temp => [1 0 1 0 1] - then add this to the current counts
%

arenas2 = arenas;

all_worm_I = sub2ind(arenaSize, sys, sxs);

temp = zeroArena;




tic
% Compute the time spent at each point for the path(s).
for i = 1:n_frames

    % Is there a skeleton for this frame?
    if isnan(sxs(1,i))
        continue;
    end

    % Compute the worm.
    wormI = sub2ind(arenaSize, sys(:,i), sxs(:,i));
    
    % Compute the time at each point for the worm points path(s).
    for j = 1:length(arenas)
        
        % Compute the unique worm points.
        wormPointsI = unique(wormI(s_points{j}));
        
        % Integrate the path.
        arenas{j}(wormPointsI) = arenas{j}(wormPointsI) + 1;
    end
end
toc

%Start of new code ...
%--------------------------------------------------------------------------
frames_run = any(all_worm_I);

tic
for iPoint = 1:n_points
   
    s_indices = s_points{iPoint};
    
    
    
    for iFrame
    
end
toc

keyboard

% Correct the y-axis (from image space).
for i = 1:length(arenas)
    arenas{i} = flipud(arenas{i});
end

% Organize the arena/path time(s).
durations(length(arenas)).indices = [];
durations(length(arenas)).times = [];
for i = 1:length(arenas)
    durations(i).indices = find(arenas{i} > 0);
    durations(i).times = double(arenas{i}(durations(i).indices)) / fps;
end

duration_struct = h__buildOutput(arena,durations);

end

function duration_struct = h__buildOutput(arena,durations)
duration_struct = struct( ...
    'arena',    arena, ...
    'worm',     durations(1), ...
    'head',     durations(2), ...
    'midbody',  durations(3), ...
    'tail',     durations(4));

end
