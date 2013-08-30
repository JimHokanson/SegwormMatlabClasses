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
%       lengths - The chain-code length at each point.
%                  
%       
%
%   See also:
%   chainCodeLength2Index
%   chainCodeLengthInterp
%   computeChainCodeLengths
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

% Are the points 2 dimensional?
if ~ismatrix(points) || (size(points, 1) ~= 2 && size(points, 2) ~= 2)
    error('circComputeChainCodeLengths:PointsNot2D', ...
        'The matrix of points must be 2 dimensional');
end

% Orient the points as a N-by-2 matrix.
isTransposed = false;
if size(points, 2) ~= 2
    points = points';
    isTransposed = true;
end

%lengths = helper__oldCode(points);

d       = @(x1,x2) sqrt(abs(x1(:,1)-x2(:,1)).^2 + abs(x1(:,2)-x2(:,2)).^2);
%NOTE: This places the wrap around length as the first element
lengths = cumsum([d(points(1,:),points(end,:)); d(points(1:end-1,:),points(2:end,:))]);

% Transpose the lengths.
if isTransposed
    lengths = lengths';
end
end

function lengths = helper__oldCode(points)

% Pre-allocate memory.
lengths = double(zeros(size(points, 1), 1));

% Pre-compute values.
sqrt2 = sqrt(2);

% Measure the difference between subsequent points.
dPoints = abs(points(1,:) - points(end,:));

% No change or we walked in a straight line.
if any(dPoints == 0)
    lengths(1) = abs(dPoints(1)) + abs(dPoints(2));
    
% We walked one point diagonally.
elseif all(dPoints == 1)
    lengths(1) = sqrt2;
    
% We walked fractionally or more than one point.
else
    lengths(1) = sqrt(sum(dPoints .^ 2));
end

% Measure the chain code length.
for i = 2:length(lengths)
    
    % Measure the difference between subsequent points.
    dPoints = abs(points(i,:) - points(i - 1,:));
    
    % No change or we walked in a straight line.
    if any(dPoints == 0)
        lengths(i) = lengths(i - 1) + abs(dPoints(1)) + abs(dPoints(2));
        
    % We walked one point diagonally.
    elseif all(dPoints == 1)
        lengths(i) = lengths(i - 1) + sqrt2;
        
    % We walked fractionally or more than one point.
    else
        lengths(i) = lengths(i - 1) + sqrt(sum(dPoints .^ 2));
    end
end

end
