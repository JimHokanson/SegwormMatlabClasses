function [distances,indices,start_index] = getLinearDistances(chain_code_lengths)
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

if transpose_output
   indices   = indices';
   distances = distances';
end