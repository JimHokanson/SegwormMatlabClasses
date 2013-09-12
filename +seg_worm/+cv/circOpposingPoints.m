function oPointsI = circOpposingPoints(pointsI, startI, endI, vLength, cc_lengths)
%circOpposingPoints  Find the equivalent point indices on the opposing side
%   of a circular vector.
%
%   pointsI = circOpposingPoints(pointsI, startI, endI, vLength, cc_lengths)
%
%   The goal is to find for each point, when the circle is divided in half,
%   points that are at the opposite mirrored side. For example:
%
%                  I
%   =>   -------------
%       s             e
%        -------------  <=
%                  O
%
%   The tricky part is that the start and end could be at any index, and
%   that the lengths may not be equal on both sides. This bases nearness
%   on the cumulative side length, NOT on physical location, i.e x-y. That
%   is accomplished in circOpposingNearestPoints.
%
%   Inputs:
%       pointsI          - the point indices to find on the opposing side
%       startI           - the index in the vector where the split, between
%                          opposing sides, starts
%       endI             - the index in the vector where the split, between
%                          opposing sides, ends
%       vLength          - the vector length
%       ??? What does this mean?
%       I guess this is passed in since the chainCodeLengths might not be
%       specified
%
%       chainCodeLengths - the chain-code length at each point;
%                          if empty, the array indices are used instead
%
%   Output:
%       pointsI - the equivalent point indices on the opposing side
%
%   See also:
%   CIRCCOMPUTECHAINCODELENGTHS
%   seg_worm.cv.circOpposingNearestPoints
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
if endI - startI < 2 || startI + vLength - endI < 2
    oPointsI = [];
    return;
end

if ~exist('cc_lengths','var') || isempty(cc_lengths)
    cc_lengths = (1:vLength)';
elseif vLength ~= length(cc_lengths);
    error('length mismatch, check code')
end


%Basically we have a set of indices, along with start and end indices
%
%   Then we have a set of lengths
%
%   Given a set of points, and the corresponds lengths, calculate
%   the index on the other side
%
%   ???? Are they flipped?
%
%   I think so, the interpretation is that they should be the same
%   amount along the worm ...

%circOpposingPoints(pointsI, startI, endI, vLength, cc_lengths)

%0 1 3 5  <- 'x' location
%5 6 8 13 <- cc_lengths 

side1_length = cc_lengths(endI) - cc_lengths(startI);
side2_length = cc_lengths(end)  - side1_length; %The remainder

on_side1_mask = pointsI > startI & pointsI < endI;
on_side2_mask = pointsI < startI | pointsI > endI;

%Convert each point to a normalized length
normalized_x1 = cc_lengths(startI:endI);
normalized_x2 = [cc_lengths(endI:end); cc_lengths(1:startI)+cc_lengths(end)];

side1_indices = startI:endI;
side2_indices = [endI:length(cc_lengths) 1:startI];

%Using the fixed indices our interpolation becomes easier
%Any point that is greater than length(cc_lengths)
%can be mapped back to the start of the indices
%i.e. length(cc_lengths)+1 => 1
side2_indices_fixed = endI:(endI+length(normalized_x2)-1);

%Actual normalization
normalized_x1 = (normalized_x1 - normalized_x1(1))./side1_length;
normalized_x2 = (normalized_x2 - normalized_x2(1))./side2_length;

%NOTE: I think #2 needs to be flipped ...
%This make it so that we mirror instead of doing a 180% match
normalized_x2 = 1-normalized_x2;

%The merger is done so that we can easily
%line up the requested points with their normalized lengths
%This is easy for side 1, but is tricky for side 2 points
%that may have moved from the front of the original series
%to the end of the normalized series ...
%
%In other words, when we extracted points from 2 for normalization
%the concatentation could have taken two separate chunks and put them
%together. Our variable pointsI only works relative to the original
%vector, not the partial vectors, so we need to merge everything back
%together again to make filling oPointsI easier.
normalized_x_merged = zeros(size(cc_lengths));
normalized_x_merged(side1_indices) = normalized_x1;
normalized_x_merged(side2_indices) = normalized_x2;

oPointsI = pointsI;

%Perform interpolation to map one side to the other
%--------------------------------------------------------------------------
%Side 1 to 2
if any(on_side1_mask)
    side1_points_x = normalized_x_merged(pointsI(on_side1_mask));

    %:/ all of this needs to be sorted, so we'll negate both sides ...
    F1 = griddedInterpolant(-1*normalized_x2,side2_indices_fixed,'linear');
    temp = round(F1(-1*side1_points_x));

    %Move any wrapped points back to the correct index 
    n_cc_lengths = length(cc_lengths);
    mask         = temp > n_cc_lengths;
    temp(mask)   = temp(mask) - n_cc_lengths;

    oPointsI(on_side1_mask) = temp;
end

%Side 2 to 1
if any(on_side2_mask)
    side2_points_x = normalized_x_merged(pointsI(on_side2_mask));
    F2 = griddedInterpolant(normalized_x1,side1_indices,'linear');
    oPointsI(on_side2_mask) = round(F2(side2_points_x));
end

% % % %DEBUGGGING
% % % %--------------------------------------------------------------------------
% % points2 = helper__oldCode(pointsI, startI, endI, vLength, cc_lengths);
% % % %--------------------------------------------------------------------------
% % % %
% % % %For debugging
% % % %----------------------------------
% % % % disp(oPointsI')
% % % % disp(points2')
% % % % 
% % if any(abs(oPointsI - points2) > 1)
% %     error('Mismatch between old code and new code')
% % end

end

function pointsI = helper__oldCode(pointsI, startI, endI, vLength, cc_lengths)



% Separate the points onto sides.
% Note: ignore start and end points, they stay the same.
% Side1 always goes from start to end in positive, index increments.
% Side2 always goes from start to end in negative, index increments.
side1  = pointsI > startI & pointsI < endI;
side2  = pointsI < startI | pointsI > endI;
side12 = side1 | side2;


%Compute the length of sides 1 & 2
%--------------------------------------------------------------------------
% Compute the size of side 1.
start_Ip1 = startI + 1;
end_Im1   = endI - 1;
%This isn't correct, should just be cc_lengths(endI) - cc_lengths(startI);
%i.e. consider a degenerate case of endI = startI + 1, then we cover
%no distance
sSize1    = cc_lengths(end_Im1) - cc_lengths(start_Ip1);

% Compute the size of side 2.
start2I = startI - 1;
if start2I < 1
    start2I = vLength;
end
end2I = endI + 1;
if end2I > vLength
    end2I = 1;
end
if start2I < end2I
    sSize2 = cc_lengths(start2I) + cc_lengths(vLength) - cc_lengths(end2I);
else % one of the ends wrapped
    sSize2 = cc_lengths(start2I) - cc_lengths(end2I);
end

% Compute the scale between sides.
scale_1to2 = sSize2 / sSize1;
scale_2to1 = sSize1 / sSize2;

% Find the distance of the side 1 points from the start, scale them for
% side 2, then find the equivalent point, at the scaled distance
% from the start, on side 2.
pointsI(side1) = cc_lengths(start2I) - (cc_lengths(pointsI(side1)) - cc_lengths(start_Ip1))*scale_1to2;

% Find the distance of the side 2 points from the start, scale them for
% side 1, then find the equivalent point, at the scaled distance
% from the start, on side 1.
minPoints2 = pointsI(side2) <= start2I;
minSide2   = false(length(pointsI),1);
maxSide2   = false(length(pointsI),1);
minSide2(side2) = minPoints2;
maxSide2(side2) = ~minPoints2;
pointsI(minSide2) = cc_lengths(start_Ip1) + ...
    (cc_lengths(start2I) - cc_lengths(pointsI(minSide2))) ...
    * scale_2to1;
pointsI(maxSide2) = cc_lengths(start_Ip1) + ...
    (cc_lengths(start2I) + cc_lengths(vLength) - ...
    cc_lengths(pointsI(maxSide2))) * scale_2to1;




% Correct any wrapped points.
wrap(side12)  = pointsI(side12) < 0;
pointsI(wrap) = pointsI(wrap) + cc_lengths(vLength);
wrap(side12)  = pointsI(side12) > cc_lengths(vLength);
pointsI(wrap) = pointsI(wrap) - cc_lengths(vLength);

% Translate the chain-code lengths to indices.
pointsI(side12) = seg_worm.cv.chainCodeLength2Index(pointsI(side12), cc_lengths);


end