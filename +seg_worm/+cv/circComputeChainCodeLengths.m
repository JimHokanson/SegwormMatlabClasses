function lengths = circComputeChainCodeLengths(points)
%circComputeChainCodeLengths Compute the chain-code length, at each point,
%   for a circularly-connected, continuous line of points.
%
%   lengths = seg_worm.cv.circComputeChainCodeLengths(points)
%
%   Input:
%       points - the circularly-connected, continuous line of points on
%                which to measure the chain-code length
%
%   Output:
%       lengths - [n x 1], The chain-code length at each point. The first element
%               is the distance from the first to the last point.
%                  
%       
%
%   See also:
%   chainCodeLength2Index
%   chainCodeLengthInterp
%   seg_worm.cv.computeChainCodeLengths

% Are the points 2 dimensional?
if ~ismatrix(points) || (size(points, 1) ~= 2 && size(points, 2) ~= 2)
    error('circComputeChainCodeLengths:PointsNot2D', ...
        'The matrix of points must be 2 dimensional');
end

% Orient the points as a N-by-2 matrix.
is_transposed = false;
if size(points, 2) ~= 2
    points = points';
    is_transposed = true;
end

d       = @(x1,x2) sqrt(abs(x1(:,1)-x2(:,1)).^2 + abs(x1(:,2)-x2(:,2)).^2);
%NOTE: This places the wrap around length as the first element
lengths = cumsum([d(points(1,:),points(end,:)); d(points(1:end-1,:),points(2:end,:))]);

% Transpose the lengths.
if is_transposed
    lengths = lengths';
end
end