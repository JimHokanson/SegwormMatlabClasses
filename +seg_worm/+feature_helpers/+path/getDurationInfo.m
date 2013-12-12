function duration_struct = getDurationInfo(sx, sy, widths, fps)
%getDurationInfo  Compute the time spent at each point along the worm path.
%
%   duration_struct = seg_worm.feature_helpers.path.getDurationInfo(sx, sy, widths, fps)
%
%   Various parts of the worm are discretized into grid locations. The
%   discretization is based on the average worm width (excludes neck and
%   hips :/ ). Once the worm skeletons have been scaled based on the width,
%   they are rounded to integer values. For each frame, all unique location
%   values increment a counter at those locations by 1. The end result is
%   that for each unit in the "arena" there is a count of how many frames
%   some part of the worm was in that location.
%
%
%   Inputs:
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
%               indices = [1  n_non_zero] the indices for the non-zero time
%                       points in the arena matrix
%               times   = [1  n_non_zero] how long some part of the worm
%                       body was at each of the indices



%Compute the scale
%--------------------------------------------------------------------------
SI = seg_worm.skeleton_indices;

% Compute the skeleton points.
s_points = {SI.ALL_INDICES SI.HEAD_INDICES SI.MID_INDICES SI.TAIL_INDICES};
n_points = length(s_points);

%??? - why scale this ????, why not just use microns?
mean_width = nanmean(mean(widths([s_points{2:4}],:),1),2);
scale = sqrt(2)/mean_width;
%NOTE: The old code omitted the widths at the neck and hips, I'm also doing
%that here, which is why we use the head, midbody, and tail indices,
%instead of just taking the mean of the widths. This distinction seems not
%functionally useful but I'll keep it for now ...

NAN_cell  = repmat({NaN},1,n_points);
durations = struct('indices',NAN_cell,'times',NAN_cell);

% The skeletons are empty.
if isempty(sx) || all(isnan(sx(:)))
    arena.height = NaN;
    arena.width = NaN;
    arena.min.x = NaN;
    arena.min.y = NaN;
    arena.max.x = NaN;
    arena.max.y = NaN;
    
    duration_struct = h__buildOutput(arena,durations);
    return;
end

% Scale the skeleton and translate so that the minimum values are at 1
%-------------------------------------------------------------------------
sxs = round(sx*scale);
sys = round(sy*scale);

xScaledMin = min(sxs(:));
xScaledMax = max(sxs(:));
yScaledMin = min(sys(:));
yScaledMax = max(sys(:));
sxs = sxs - xScaledMin + 1;
sys = sys - yScaledMin + 1;


% Construct the empty arena(s).
arenaSize = [yScaledMax - yScaledMin + 1, xScaledMax - xScaledMin + 1];

%Organize the arena size.
%------------------------------------------------
arena.height = arenaSize(1);
arena.width  = arenaSize(2);
arena.min.x  = min(sx(:));
arena.min.y  = min(sy(:));
arena.max.x  = max(sx(:));
arena.max.y  = max(sy(:));
%--------------------------------------------------------------------------

arenas = h__populateArenas(arenaSize, sys, sxs, s_points);


% Organize the arena/path time(s).
for iPoint = 1:n_points
    non_empty_arena_indices   = find(arenas{iPoint} > 0);
    durations(iPoint).indices = non_empty_arena_indices;
    durations(iPoint).times   = arenas{iPoint}(non_empty_arena_indices) / fps;
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

function arenas = h__populateArenas(arena_size, sys, sxs, s_points)
%
%
%   
%
%

%NOTE: All skeleton points have been rounded to integer values for
%assignment to the matrix based on their values being treated as indices

%Convert to linear indices for assignment.
all_worm_I   = sub2ind(arena_size, sys, sxs);
frames_run   = find(any(all_worm_I)); %NOTE: any(NaN) is false, all(NaN) is true
n_frames_run = length(frames_run);

%1 area for each set of skeleton indices
arenas    = cell(n_points,1);

for iPoint = 1:n_points
   
    temp_arena = zeros(arena_size);
    s_indices  = s_points{iPoint};
    
    for iFrame = 1:n_frames_run
       cur_frame   = frames_run(iFrame);
       cur_indices = all_worm_I(s_indices,cur_frame);
       

       %Approach: we only want to increment 1 for each unique value, but
       %assuming that the right hand side is done before any assigments are
       %made, redundant assignments result in only having each unique value
       %incremented by 1. Previously this code used unique() which is much
       %slower.
       %
       %i.e. a = [0 0 0 0 0]
       %     b = [1 3 1 3 5] %NOTE: We have two each of 1 & 3
       %
       %     a(b) = a(b) + 1 => [1 0 1 0 1]
       %
       %    I assume the computer does the assignment:
       %    a(1) = 1, twice 
       %
       %    and not:
       %    a(1) = a(1) + 1 twice
       %    
       %    same for 3:
       %    a(3) = 1 , NOT a(3) = a(3) + 1
       %
       %    This allows us to avoid computing the unique set of indices
       %    before doing the calculation:
       %    i.e., we avoid b = unique(b)
       %
       temp_arena(cur_indices) = temp_arena(cur_indices) + 1;
    end
    
    % Correct the y-axis (from image space).
    %???? Hold over from old code
    arenas{iPoint} = temp_arena(end:-1:1,:);
end

end