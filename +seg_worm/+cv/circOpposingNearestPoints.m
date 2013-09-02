function pointsI = circOpposingNearestPoints(points_I_in, x, startI, endI, ...
    searchLength, cc_lengths)
%circOpposingNearestPoints Find the nearest equivalent point indices on the
%   opposing side (within a search window) of a circular vector.
%
%   pointsI = circOpposingNearestPoints(...
%                   pointsI, X, startI, endI, searchLength, *cc_lengths)
%
%   ???? - it is not clear how this would work close to the head and tail
%   ...
%
%
%   Inputs:
%       pointsI          - the point indices to find on the opposing side
%       x                - [n x 1] the circularly connected vector on which the
%                          points lie
%       startI           - the index in the vector where the split, between
%                          opposing sides, starts
%       endI             - the index in the vector where the split, between
%                          opposing sides, ends
%       searchLength     - the search length, on either side of a directly
%                          opposing point, to search for the nearest point
%       chainCodeLengths - the chain-code length at each point;
%                          if empty, the array indices are used instead
%
%   Output:
%       pointsI - the equivalent point indices on the opposing side
%
%   See also:
%   CIRCOPPOSINGPOINTS
%   CIRCNEARESTPOINTS
%   CIRCCOMPUTECHAINCODELENGTHS
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software;
% you must reproduce all copyright notices and other proprietary
% notices on any copies of the Software.


% Re-order the start and end to make life simple.
if startI > endI
    [startI,endI] = deal(endI,startI);
end

% The points are degenerate.
if endI - startI < 2 || startI + size(x,1) - endI < 2
    pointsI = [];
    return;
end

pointsI = points_I_in;

% Are there chain-code lengths?
if ~exist('cc_lengths','var') || isempty(cc_lengths)
    cc_lengths = (1:size(x,1))';
end

% Compute the opposing points.
oPointsI = seg_worm.cv.circOpposingPoints(pointsI, startI, endI, size(x,1), cc_lengths);

%This is temporary, for debugging purposes
final_old = helper_oldCode(x,pointsI,oPointsI,startI,endI,searchLength, cc_lengths);

%The goal is for each opposite index to get a range of distances to search
%so that we can get the opposite location that is closest based on x-y, not
%based on the cumulative length of the side (what the above function does)


on_a_side_mask = oPointsI ~= startI & oPointsI ~= endI;
oPointsI(~on_a_side_mask) = [];

n_side_points = length(oPointsI);

minOPointsI = NaN(n_side_points,1);
maxOPointsI = NaN(n_side_points,1);

[distances,indices,~,x_locations] = seg_worm.util.getLinearDistances(cc_lengths);

side1_mask = oPointsI > startI & oPointsI < endI;
side2_mask = oPointsI < startI | oPointsI > endI;

if any(side1_mask)
    target_1_indices = indices;
    %We apply a fix so that if we are outside the valid
    %index range on either side, we get bumped to the boundary
    %
    %??? - change for equal as well????
    %I think so, since normally we bump to extremes
    %
    %   2 2 2 s 1 1 1 1 1 e 2 2 2
    %     L         R
    %   We move the left bound to:
    %           L
    %
    target_1_indices(indices >= endI)   = endI-1;
    target_1_indices(indices <= startI) = startI + 1;
    
    [minOPointsI(side1_mask),maxOPointsI(side1_mask)] = ...
        helper__updateMinMax(distances,x_locations,target_1_indices,oPointsI,side1_mask,searchLength);
    
end

if any(side2_mask)    
    if startI == 1
        right_value_2 = length(cc_lengths);
    else
        right_value_2 = startI-1;
    end
    
    if endI == length(cc_lengths)
        left_value_2 = 1;
    else
        left_value_2 = endI + 1;
    end
        
    %Anything that is a side 1 index, move to the bounds
    %This one is a bit tricky to know which index is closer to our index
    %
    %   This bit of code determines that, and we assign
    %   the edge values, right_value_2 and left_value_2 accordingly
    x_start = x_locations(startI);
    x_end   = x_locations(endI);
    x_half  = x_start + 0.5*(x_end - x_start);
    [~,half_I] = min(abs(x_locations - x_half));
    
    side1_mask_or_SE = indices >= startI & indices < endI;
    
    target_2_indices = indices;
    target_2_indices(side1_mask_or_SE & indices >= half_I) = right_value_2;
    target_2_indices(side1_mask_or_SE & indices < half_I)  = left_value_2;
    
    [minOPointsI(side2_mask),maxOPointsI(side2_mask)] = ...
        helper__updateMinMax(distances,x_locations,target_2_indices,oPointsI,side2_mask,searchLength);
    
    
end

pointsI(on_a_side_mask) = seg_worm.cv.circNearestPoints(x(pointsI(on_a_side_mask),:), minOPointsI, maxOPointsI, x);

if ~isequal(pointsI,final_old)
   error('Mismatch between my implementation and old code version') 
end

end

function [minO,maxO] = helper__updateMinMax(distances,x_locations,target_indices,oPointsI,mask,searchLength)

one_to_N = (1:length(distances))';

x_locations_1      = x_locations(oPointsI(mask));

left_dist_targets  = x_locations_1 - searchLength;
right_dist_targets = x_locations_1 + searchLength;

minO = target_indices(ceil(interp1(distances,one_to_N,left_dist_targets)));
maxO = target_indices(floor(interp1(distances,one_to_N,right_dist_targets)));
end


function pointsI = helper_oldCode(x,pointsI,oPointsI,startI,endI,searchLength, cc_lengths)

import seg_worm.cv.*


% Separate the points onto sides.
% Note: ignore start and end points, they stay the same.
% Side1 always goes from start to end in positive, index increments.
% Side2 always goes from start to end in negative, index increments.
side12 = oPointsI ~= startI & oPointsI ~= endI;
oPointsI(~side12) = [];
side1 = oPointsI > startI & oPointsI < endI;
side2 = oPointsI < startI | oPointsI > endI;

% Compute the start indices.
% Note: we checked for degeneracy; therefore, only one index can wrap.
is2Wrap = false;
start1  = startI + 1;
start2  = startI - 1;
if start2 < 1
    start2  = start2 + size(x,1);
    is2Wrap = true;
end

% Compute the end indices.
end1 = endI - 1;
end2 = endI + 1;
if end2 > size(x,1)
    end2    = end2 - size(x,1);
    is2Wrap = true;
end

% Compute the minimum search points on side 2 (for the search intervals opposite side 1).
minOPointsI(side1) = cc_lengths(oPointsI(side1)) - searchLength;
wrap               = false(size(side1));
wrap(side1)        = minOPointsI(side1) < cc_lengths(1);
minOPointsI(wrap)  = minOPointsI(wrap) + cc_lengths(end);
wrap               = false(size(side1));

wrap(side1) = minOPointsI(side1) < cc_lengths(start1) | ...
    minOPointsI(side1) > cc_lengths(end1);

minOPointsI(wrap)    = start1;
notWrap              = side1 & ~wrap;
minOPointsI(notWrap) = chainCodeLength2Index(minOPointsI(notWrap),cc_lengths);

% Compute the maximum search points on side 2 (for the search intervals opposite side 1).
maxOPointsI(side1) = cc_lengths(oPointsI(side1)) + searchLength;
wrap               = false(size(side1));
wrap(side1)        = maxOPointsI(side1) > cc_lengths(end);
maxOPointsI(wrap)  = maxOPointsI(wrap) - cc_lengths(end);
wrap               = false(size(side1));

wrap(side1) = maxOPointsI(side1) < cc_lengths(start1) | maxOPointsI(side1) > cc_lengths(end1);

maxOPointsI(wrap)    = end1;
notWrap              = side1 & ~wrap;
maxOPointsI(notWrap) = chainCodeLength2Index(maxOPointsI(notWrap),cc_lengths);

% Compute the minimum search points on side 1 (for the search intervals opposite side 2).
minOPointsI(side2) = cc_lengths(oPointsI(side2)) - searchLength;
wrap               = false(size(side2));
wrap(side2)        = minOPointsI(side2) < cc_lengths(1);
minOPointsI(wrap)  = minOPointsI(wrap) + cc_lengths(end);
wrap               = false(size(side2));
if is2Wrap
    wrap(side2) = minOPointsI(side2) > cc_lengths(start2) | minOPointsI(side2) < cc_lengths(end2);
else
    wrap(side2) = minOPointsI(side2) > cc_lengths(start2) & minOPointsI(side2) < cc_lengths(end2);
end
minOPointsI(wrap)    = end2;
notWrap              = side2 & ~wrap;
minOPointsI(notWrap) = chainCodeLength2Index(minOPointsI(notWrap), cc_lengths);

% Compute the maximum search points on side 1 (for the search intervals opposite side 2).
maxOPointsI(side2) = cc_lengths(oPointsI(side2)) + searchLength;
wrap               = false(size(side2));
wrap(side2)        = maxOPointsI(side2) > cc_lengths(end);
maxOPointsI(wrap)  = maxOPointsI(wrap) - cc_lengths(end);
wrap               = false(size(side2));

if is2Wrap
    wrap(side2) = maxOPointsI(side2) > cc_lengths(start2) | ...
        maxOPointsI(side2) < cc_lengths(end2);
else
    wrap(side2) = maxOPointsI(side2) > cc_lengths(start2) & ...
        maxOPointsI(side2) < cc_lengths(end2);
end

maxOPointsI(wrap)    = start2;
notWrap              = side2 & ~wrap;
maxOPointsI(notWrap) = chainCodeLength2Index(maxOPointsI(notWrap),cc_lengths);
 
% disp(minOPointsI)
% disp(maxOPointsI)

% Search for the nearest points.
pointsI(side12) = circNearestPoints(x(pointsI(side12),:), minOPointsI, maxOPointsI, x);


end
