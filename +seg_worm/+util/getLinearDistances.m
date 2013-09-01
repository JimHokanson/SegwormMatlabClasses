function [distances,indices,start_index,x_locations] = getLinearDistances(chain_code_lengths)
%
%   This is a helper function which returns distances such that neighbor
%   operations are possible in both directions without worrying about edge
%   effects.
%   
%
%   chain_code_lengths
%   1st element - circular distance from first to last point
%   2nd element - distance from 1st to 2nd point + value at index 1
%
%
%   points      = [0 1 3 6 10 15]
%   point_diffs = [15 1 2 3 4  5]
%   cc_length   = [15 16 18 21 25 30]
%
%   See Also:
%   seg_worm.cv.circCurvatures
%   seg_worm.util.peaksCircDist

if iscolumn(chain_code_lengths)
    chain_code_lengths = chain_code_lengths';
    transpose_output = true;
else
    transpose_output = false;
end

%[-30,-29,-27,-24,-20,-15]
reverse_distance = fliplr(cumsum([-1*chain_code_lengths(1) diff(chain_code_lengths(end:-1:1))]));

%[0,1,3,6,10,15]
normal_data = [0 chain_code_lengths(2:end)-chain_code_lengths(1)];

%[30,31,33,36,40,45]
extended_data = normal_data+chain_code_lengths(end);

%[-30,-29,-27,-24,-20,-15,0,1,3,6,10,15,30,31,33,36,40,45]
distances = [reverse_distance normal_data extended_data];

n_points = length(chain_code_lengths);

indices = repmat(1:n_points,1,3);

start_index = n_points + 1;

x_locations = normal_data;

if transpose_output
   indices   = indices';
   distances = distances';
   x_locations = x_locations';
end

end

%   A bit more on the algorithm: 
%   Every point is given an absolute location based on distance. We then
%   replicate the points once on the left and the right, to account for the
%   edges. For each new peak, we turn off all points that are within a
%   certain distance of that peak.
%
%   Worked example (see below). Let's say we choose c, which is at 4. Then
%   let's say we do +/- 5, this puts our edges at -1 and 9. From the chart
%   below, we see that this encompasses a through c, so a&b can not be used
%   as peaks since they are too close to c. Note, that if a or b had
%   already been chosen, then c would have been ineleligble as a peak. If
%   we expand the max distance to being +/- 6, instead of 5, then c would
%   just be touching e and d. d is obvious but e would be touched because
%   it is at -2 (4 - 6 => -2)
%
%Consider the following
%               1 2 3   4    5  index
%------------------------------------------------------------------
%               a b c   d   e   'a' <= points
%                1 3 6   7    2     <= distance between cells
%                1 4 10  17   19    <= chain_code_lengths
%               0 1 4   10  19      <= absolute location 
%
%             x
%   a   b   c   d   e   a   b   c   d   e   a   b   c   d   e
%     1   3   6   7   2   1   3   6   7   2   1   3   6   7   2       
%      -18  -15 -9  -2  0   1   4   10  17  19  20  23  29  36     
%
%       2   3   4   5   1   2   3   4   5   1   2   3   4   5

%These few lines are a mess, but make more sense if you look at the example above
%--------------------------------------------------------------------