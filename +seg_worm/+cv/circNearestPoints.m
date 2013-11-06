function nearI = circNearestPoints(points, minI, maxI, x)
%circNearestPoints  For each point, find the nearest corresponding point
%   within an interval of circularly-connected search points (based on
%   distance - i.e. x,y location)
%
%   nearI = circNearestPoints(points, minI, maxI, x)
%
%   This looks for the closest x-y point on the opposite
%   side within a given range.
%
%
%
%   Inputs:
%       points - [n x 2] the point coordinates from which the distance is measured
%       minI   - [n x 2] the minimum indices of the intervals
%       maxI   - [n x 2] the maximum indices of the intervals
%       x      - the circularly-connected, point coordinates on which the
%                search intervals lie
%
%   Output:
%       nearI - the indices of the nearest points
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

% Are the points 2 dimensional?
if ~ismatrix(points) || (size(points, 1) ~= 2 && size(points, 2) ~= 2)
    error('circNearestPoints:PointsNot2D', ...
        'The matrix of points must be 2 dimensional');
end

% Are the search points 2 dimensional?
if length(size(x)) ~=2 || (size(x, 2) ~= 2 && size(x, 2) ~= 2)
    error('circNearestPoints:XNot2D', ...
        'The circularly-connected search points must be 2 dimensional');
end

% Orient the points as a N-by-2 matrix.
isTransposed = false;
if size(points, 2) ~= 2
    points = points';
    isTransposed = true;
end

% Orient the search points as a N-by-2 matrix.
if size(x, 2) ~= 2
    x = x';
end

n_points = size(points,1);
nearI    = NaN(n_points,1);

% Search for the nearest points.
for iPoint = 1:n_points
    x1 = points(iPoint,:);
    
    if minI(iPoint) <= maxI(iPoint)
        % The interval is continuous.
        x2 = x(minI(iPoint):maxI(iPoint),:);
        [~, temp] = min(pdist2(x1,x2));
        
        nearI(iPoint) = temp + minI(iPoint) - 1;
    else
        % The interval wraps.
        x2 = x(minI(iPoint):end,:);
        x3 = x(1:maxI(iPoint),:);
        
        [mag1, nearI1] = min(pdist2(x1,x2));
        [mag2, nearI2] = min(pdist2(x1,x3));
        
        % Which point is nearest?
        if mag1 <= mag2
            nearI(iPoint) = nearI1 + minI(iPoint) - 1;
        else
            nearI(iPoint) = nearI2;
        end
    end
end

% Transpose the point indices.
if isTransposed
    nearI = nearI';
end
end
