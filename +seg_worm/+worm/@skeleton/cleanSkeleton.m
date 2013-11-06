function [cSkeleton,cWidths] = cleanSkeleton(obj,skeleton, widths, wormSegSize)
%CLEANSKELETON Clean an 8-connected skeleton by removing any overlap and
%interpolating any missing points.
%
%   [CSKELETON] = CLEANSKELETON(SKELETON)
%
%   Note: the worm's skeleton is still rough. Therefore, index lengths, as
%         opposed to chain-code lengths, are used as the distance metric
%         over the worm's skeleton.
%
%   Input:
%       skeleton    - the 8-connected skeleton to clean
%       widths      - the worm's contour widths at each skeleton point
%       wormSegSize - the size (in contour points) of a worm segment.
%                     Note: the worm is roughly divided into 24 segments
%                     of musculature (i.e., hinges that represent degrees
%                     of freedom) on each side. Therefore, 48 segments
%                     around a 2-D contour.
%                     Note 2: "In C. elegans the 95 rhomboid-shaped body
%                     wall muscle cells are arranged as staggered pairs in
%                     four longitudinal bundles located in four quadrants.
%                     Three of these bundles (DL, DR, VR) contain 24 cells
%                     each, whereas VL bundle contains 23 cells." -
%                     www.wormatlas.org
%
%   Output:
%       cSkeleton - the cleaned skeleton (no overlap & no missing points)
%       cWidths   - the cleaned contour widths at each skeleton point
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.


% If a worm touches itself, the cuticle prevents the worm from folding and
% touching adjacent pairs of muscle segments; therefore, the distance
% between touching segments must be, at least, the length of 2 muscle
% segments.
maxSkeletonOverlap = 2 * wormSegSize;

%Removal of small loops
%--------------------------------------------------------------------------
%[skeleton2,widths2] = helper__removeSmallLoopsOld(skeleton,widths,maxSkeletonOverlap);

[skeleton,widths] = helper__removeSmallLoops(skeleton,widths,maxSkeletonOverlap);

% if ~isequal(skeleton,skeleton2) || ~isequal(widths,widths2)
%    error('Mismatch between old and new')
% end

%Headling of missing points
%--------------------------------------------------------------------------
%[cSkeleton2,cWidths2] = helper__healMissingPointsOld(skeleton,widths);

[cSkeleton,cWidths] = helper__healMissingPoints(skeleton,widths);

% if ~isequal(cSkeleton,cSkeleton2) || ~isequal(cWidths,cWidths2)
%    error('Mismatch between old and new')
% end


%Anti-aliasing the skeleton
%--------------------------------------------------------------------------
%[cSkeleton2,cWidths2] = helper__antiAliasOld(cSkeleton,cWidths);

[cSkeleton,cWidths]   = helper__antiAlias(cSkeleton,cWidths);
%toc

% if ~isequal(cSkeleton,cSkeleton2) || ~isequal(cWidths,cWidths2)
%    error('Mismatch between old and new')
% end

end

function [cSkeleton,cWidths] = helper__antiAliasOld(cSkeleton,cWidths)
% Anti alias.
%--------------------------------------------------------------------------
keep = 1:size(cSkeleton, 1); % points to keep
i    = 1;
endI = size(cSkeleton, 1) - 1;
while i < endI
    
    % Smooth any stairs.
    nextI = i + 2;
    if abs(cSkeleton(i,1) - cSkeleton(nextI,1)) <= 1 && abs(cSkeleton(i,2) - cSkeleton(nextI,2)) <= 1
        keep(i + 1) = nan;
        
        % Advance.
        i = nextI;
        
        % Advance.
    else
        i = i + 1;
    end
end
cSkeleton = cSkeleton(~isnan(keep),:);
cWidths = cWidths(~isnan(keep));
end

function [cSkeleton,cWidths] = helper__antiAlias(cSkeleton,cWidths)
%??? - why is this not like the other anti-alias
%which compares the next 2 points ??????
%
%
%   See seg_worm.contour.cleanWorm/helper__removeUnecessaryContourPoints

% Anti alias.
%--------------------------------------------------------------------------
n_points  = size(cSkeleton,1);
keep_mask = true(1,n_points);

i = 1;

while i < n_points - 1
    
    % Smooth any stairs.
    nextI = i + 2;
    if abs(cSkeleton(i,1) - cSkeleton(nextI,1)) <= 1 && abs(cSkeleton(i,2) - cSkeleton(nextI,2)) <= 1
        keep_mask(i + 1) = false;
        
        % Advance.
        i = nextI;
    else
        % Advance.
        i = i + 1;
    end
end

%NOTE: The above code doesn't compare the last two points. It is possible
%that the last two points are the same, which causes problems ...
if cSkeleton(end-1,1) == cSkeleton(end,1) && cSkeleton(end-1,2) == cSkeleton(end,2)
   keep_mask(end-1) = false; 
end

if cSkeleton(end-2,1) == cSkeleton(end,1) && cSkeleton(end-2,2) == cSkeleton(end,2)
   keep_mask(end-2) = false;  
end


cSkeleton = cSkeleton(keep_mask,:);
cWidths   = cWidths(keep_mask);
end



function [skeleton,widths] = helper__removeSmallLoopsOld(skeleton,widths,maxSkeletonOverlap)

% Remove small loops.
keep = 1:size(skeleton, 1); % points to keep
[~, pSort] = sortrows(skeleton); % the sorted points
[~, iSort] = sort(pSort); % index -> sorted point index
s1I = 1; % the first index for the skeleton loop
while s1I < length(pSort)
    
    % Find small loops.
    % Note: distal, looped sections are most likely touching;
    % therefore, we don't remove these.
    if ~isnan(keep(s1I))
        minI = s1I; % the minimum index for the loop
        maxI = s1I; % the maximum index for the loop
        
        % Search backwards.
        if iSort(s1I) > 1
            pI = iSort(s1I) - 1; % the index for the sorted points
            s2I = pSort(pI); % the second index for the skeleton loop
            dSkeleton = abs(skeleton(s1I,:) - skeleton(s2I,:));
            while any(dSkeleton <= 1)
                if s2I > s1I && ~isnan(keep(s2I)) && ...
                        all(dSkeleton <= 1) && ...
                        abs(s1I - s2I) < maxSkeletonOverlap
                    minI = min(minI, s2I);
                    maxI = max(maxI, s2I);
                end
                
                % Advance the second index for the skeleton loop.
                pI = pI - 1;
                if pI < 1
                    break;
                end
                s2I = pSort(pI);
                dSkeleton = abs(skeleton(s1I,:) - skeleton(s2I,:));
            end
        end
        
        % Search forwards.
        if  iSort(s1I) < length(pSort)
            pI = iSort(s1I) + 1; % the index for the sorted points
            s2I = pSort(pI); % the second index for the skeleton loop
            dSkeleton = abs(skeleton(s1I,:) - skeleton(s2I,:));
            while any(dSkeleton <= 1)
                if s2I > s1I && ~isnan(keep(s2I)) && ...
                        all(dSkeleton <= 1) && ...
                        abs(s1I - s2I) < maxSkeletonOverlap
                    minI = min(minI, s2I);
                    maxI = max(maxI, s2I);
                end
                
                % Advance the second index for the skeleton loop.
                pI = pI + 1;
                if pI > length(pSort)
                    break;
                end
                s2I = pSort(pI);
                dSkeleton = abs(skeleton(s1I,:) - skeleton(s2I,:));
            end
        end
        
        % Remove small loops.
        if minI < maxI
            
            % Remove the overlap.
            if isequal(skeleton(minI,:), skeleton(maxI,:))
                keep((minI + 1):maxI) = nan;
                widths(minI) = min(widths(minI:maxI));
                
            % Remove the loop.
            elseif minI < maxI - 1
                keep((minI + 1):(maxI - 1)) = nan;
                widths(minI) = min(widths(minI:(maxI - 1)));
                widths(maxI) = min(widths((minI + 1):(maxI)));
            end
            
        end
        
        % Advance the first index for the skeleton loop.
        if s1I < maxI
            s1I = maxI;
        else
            s1I = s1I + 1;
        end
        
    % Advance the first index for the skeleton loop.
    else
        s1I = s1I + 1;
    end
end
skeleton = skeleton(~isnan(keep),:);
widths = widths(~isnan(keep));

end

function [skeleton,widths] = helper__removeSmallLoops(skeleton,widths,maxSkeletonOverlap)

max_distance = sqrt(2); %NOTE: This algorithm does <= radius
IDX  = rangesearch(skeleton,skeleton,max_distance); %Stats toolbox call

n_points  = size(skeleton,1);
keep_mask = true(1,n_points);

%Update min or max if:
%1) not taken
%2) It is greater than our current point 
%   - note, this means the min should never be updated ...
%3) Distance is within maxSkeletonOverlap
%
%
%   When updating, 
%   -- if equal, then remove all points between current
%   and the max, update widths using min to max
%   -- if the difference is greater than 1, remove the points in between,
%       update widths using min to max - 1

i1 = 1;
while i1 < n_points
    if keep_mask(i1)
        cur_indices = IDX{i1};
        %Grab max of all the indices we are keeping or the same index
        %NOTE: the second max just takes care of the output not being empty
        %i.e. we have max(i1,[]) -> just give us i1
        %
        %               This gets the max of all points that are close
        i2 = max(i1,   max(cur_indices(keep_mask(cur_indices)))    );
        
        %If our max is too far away, indicating a really big loop
        %then we need to do a slower approach where we also ensure
        %that the distance is good ...
        if i2 - i1 > maxSkeletonOverlap
            %Then we need an additional check to find the point that is closest
            %Hopefully we won't normally need to do this, which is why
            %we have an extra check
            
            temp_indices = cur_indices(keep_mask(cur_indices));
            i2 = max(i1,    max(temp_indices(temp_indices - i1 < maxSkeletonOverlap))    );
        end
        
        if i2 > i1
            if isequal(skeleton(i1,:), skeleton(i2,:))
                %Remove the loop since in a very small space we get
                %back to the same location.
                keep_mask((i1 + 1):i2) = false;
                widths(i1) = min(widths(i1:i2));
            elseif i2 > i1 + 1
                % Remove the loop.
                keep_mask((i1 + 1):(i2 - 1)) = false;
                widths(i1) = min(widths(i1:(i2 - 1)));
                widths(i2) = min(widths((i1 + 1):i2));
            end
            i1 = i2;
        else
            i1 = i1 + 1;
        end
    else
        i1 = i1 + 1;
    end
end

skeleton = skeleton(keep_mask,:);
widths   = widths(keep_mask);
end

function [cSkeleton,cWidths] = helper__healMissingPointsOld(skeleton,widths)

% The head and tail have no width.
widths(1) = 0;
widths(end) = 0;

% Heal the skeleton by interpolating missing points.
cSkeleton = zeros(2 * size(skeleton, 1), 2); % pre-allocate memory
cWidths   = zeros(2 * size(skeleton, 1), 1); % pre-allocate memory
j = 1;
for i = 1:(length(skeleton) - 1)
    
    % Initialize the point differences.
    y = abs(skeleton(i + 1,1) - skeleton(i,1));
    x = abs(skeleton(i + 1,2) - skeleton(i,2));
    
    % Add the point.
    if (y == 0 || y == 1) && (x == 0 || x == 1)
        cSkeleton(j,:) = skeleton(i,:);
        cWidths(j) = widths(i);
        j = j + 1;
        
    % Interpolate the missing points.
    else
        points = max(y, x);
        y1 = skeleton(i,1);
        y2 = skeleton(i + 1,1);
        x1 = skeleton(i,2);
        x2 = skeleton(i + 1,2);
        cSkeleton(j:(j + points),1) = round(linspace(y1, y2, points + 1));
        cSkeleton(j:(j + points),2) = round(linspace(x1, x2, points + 1));
        cWidths(j:(j + points)) = round(linspace(widths(i), ...
            widths(i + 1), points + 1));
        j = j + points;
    end
end

% Add the last point.
if (cSkeleton(1,1) ~= skeleton(end,1)) || ...
        (cSkeleton(1,2) ~= skeleton(end,2))
    cSkeleton(j,:) = skeleton(end,:);
    cWidths(j) = widths(end);
    j = j + 1;
end

% Collapse any extra memory.
cSkeleton(j:end,:) = [];
cWidths(j:end) = [];

end

function [cSkeleton,cWidths] = helper__healMissingPoints(skeleton,widths)
%
%   This code looks for gaps that are larger than sqrt(2) and fills
%   them in using linear interpolation.
%   
%   NOTE: This code could be optimized a bit more but it probably isn't
%   worth it ...

% The head and tail have no width.
widths(1)   = 0;
widths(end) = 0;

% Heal the skeleton by interpolating missing points.
cSkeleton = zeros(2 * size(skeleton, 1), 2); % pre-allocate memory
cWidths   = zeros(2 * size(skeleton, 1), 1); % pre-allocate memory
j = 1;

d      = sqrt(sum(diff(skeleton,[],1).^2,2));
no_fix = d <= sqrt(2);

for i = 1:(length(skeleton) - 1)
    
    % Add the point.
    if no_fix(i)
        cSkeleton(j,:) = skeleton(i,:);
        cWidths(j)     = widths(i);
        j = j + 1;
        
    % Interpolate the missing points.
    else
        dy = abs(skeleton(i + 1,1) - skeleton(i,1));
        dx = abs(skeleton(i + 1,2) - skeleton(i,2));
        points = max(dy, dx);
        y1 = skeleton(i,1);
        y2 = skeleton(i + 1,1);
        x1 = skeleton(i,2);
        x2 = skeleton(i + 1,2);
        cSkeleton(j:(j + points),1) = round(linspace(y1, y2, points + 1));
        cSkeleton(j:(j + points),2) = round(linspace(x1, x2, points + 1));
        cWidths(j:(j + points))     = round(linspace(widths(i), widths(i + 1), points + 1));
        j = j + points;
    end
end

% Add the last point.
if (cSkeleton(1,1) ~= skeleton(end,1)) || (cSkeleton(1,2) ~= skeleton(end,2))
    cSkeleton(j,:) = skeleton(end,:);
    cWidths(j)     = widths(end);
    j = j + 1;
end

% Collapse any extra memory.
cSkeleton(j:end,:) = [];
cWidths(j:end) = [];

end