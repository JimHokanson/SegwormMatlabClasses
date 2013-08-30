function [peaks,indices] = peaksCircDist(x, dist,use_max,value_cutoff,chain_code_lengths)
%MAXPEAKSCIRCDIST Find the maximum peaks in a circular vector. The peaks
%are separated by, at least, the given distance.
%
%   [peaks,indices] = seg_worm.util.peaksCircDist(x, dist, *chain_code_lengths)
%
%   Algorithm: ???? - grab best peaks or sequential peaks start at index 1?
%
%       I think this algorithm grabs the best peaks, starting with the
%       first, then grabbing the 2nd (if we can), then the 3rd, etc ...
%
%   Inputs:
%       x                - the vector of values
%       dist             - the minimum distance between peaks
%       chain_code_lengths - the chain-code length at each index;
%                          if empty, the array indices are used instead
%
%   Outputs:
%       peaks   - extreme values
%       indices - the indices of these extreme values, ordered by
%           original index, not by value, i.e. indices(1) < indices(2)
%           always
%
%   See also:
%   CIRCCOMPUTECHAINCODELENGTHS

if ~exist('chain_code_lengths','var') || isempty(chain_code_lengths)
    % Use the array indices for length.
    chain_code_lengths = 1:length(x);
else
    keyboard
end

if use_max
    I1 = find(x > value_cutoff);
    [~,I2] = sort(-1*x(I1));
    I = I1(I2);
else
    I1     = find(x < value_cutoff);
    [~,I2] = sort(x(I1));
    I = I1(I2);
end

n_sort = length(I);
n_x = length(x);

if size(x,1) > 1
    x = x';
end


%   Algorithm: 
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
rev_data = cumsum(diff(chain_code_lengths(end:-1:1)));
end_data = chain_code_lengths(end) + chain_code_lengths(1:end-1);
relative_distances  = [rev_data(end:-1:1) 0 chain_code_lengths end_data];
offset_index_values = [2:n_x 1:n_x 1:n_x];


x_location = [0 chain_code_lengths(1:end-1)];

x_cutoff_left  = x_location - dist;
x_cutoff_right = x_location + dist; 

%Remove need for -1
%TODO: Move computeEdgeIndices to STD_LIB
%start <= distance < end

[start_I,end_I] = computeEdgeIndices(relative_distances,x_cutoff_left,x_cutoff_right);
start_I = start_I - 1;

% taken_start_I       = offset_index_values(start_I);
% taken_end_I         = offset_index_values(end_I);

is_peak_mask = false(1,n_x); 
taken_mask   = false(1,n_x); %Indicates that the index is close to or is 
%a peak and thus can not be used as a peak
for iElem = 1:n_sort
    cur_index = I(iElem);
    if ~taken_mask(cur_index)
       is_peak_mask(cur_index)  = true;
       temp_indices             = offset_index_values(start_I(cur_index):end_I(cur_index));
       taken_mask(temp_indices) = true;
    end
end

indices = find(is_peak_mask);
peaks   = x(indices);

%[peaks,indices] = helper__oldCode(x, dist, chain_code_lengths,win_size)



end

function [peaks,indices] = helper__oldCode(x, dist, chain_code_lengths,winSize) %#ok<DEFNU>
% Initialize the peaks and indices.
wins    = ceil(length(x) / winSize);
peaks   = zeros(wins, 1); % pre-allocate memory
indices = zeros(wins, 1); % pre-allocate memory

% Search for peaks.
im = 0; % the last maxima index
ie = 0; % the end index for the last maxima's search window
ip = 1; % the current, potential, max peak index
p  = x(ip); % the current, potential, max peak value
i  = 2; % the vector index
j  = 1; % the recorded, maximal peaks index
while i <= length(x)
    
    % Found a potential peak.
    if x(i) > p
        ip = i;
        p = x(i);
    end
    
    % Test the potential peak.
    if chain_code_lengths(i) - chain_code_lengths(ip) >= dist || i == length(x)
        
        % Check the untested values next to the previous maxima.
        if im > 0 && chain_code_lengths(ip) - chain_code_lengths(im) <= 2 * dist
            
            % Check the untested values next to the previous maxima.
            isMax = true;
            k = ie;
            while isMax && k > 0 && ...
                    chain_code_lengths(ip) - chain_code_lengths(k) < dist
                
                % Is the previous peak larger?
                if x(ip) <= x(k)
                    isMax = false;
                end
                
                % Advance.
                k = k - 1;
            end
            
            % Record the peak.
            if isMax
                indices(j) = ip;
                peaks(j) = p;
                j = j + 1;
            end
            
            % Record the maxima.
            im = ip;
            ie = i;
            ip = i;
            p = x(ip);
            
            % Record the peak.
        else
            indices(j) = ip;
            peaks(j) = p;
            j = j + 1;
            im = ip;
            ie = i;
            ip = i;
            p = x(ip);
        end
    end
    
    % Advance.
    i = i + 1;
end

% Collapse any extra memory.
indices(j:end) = [];
peaks(j:end) = [];

% If we have two or more peaks, we have to check the start and end for mistakes.
if j > 2
    
    % If the peaks at the start and end are too close, keep the largest or
    % the earliest one.
    if (chain_code_lengths(indices(1)) + chain_code_lengths(end) - ...
            chain_code_lengths(indices(end))) < dist
        if peaks(1) <= peaks(end)
            indices(1) = [];
            peaks(1) = [];
        else
            indices(end) = [];
            peaks(end) = [];
        end
        
        % Otherwise, check any peaks that are too close to the start and end.
    else
        
        % If we have a peak at the start, check the wrapping portion just
        % before the end.
        k = length(x);
        while chain_code_lengths(indices(1)) + ...
                chain_code_lengths(length(x)) - chain_code_lengths(k) < dist
            
            % Remove the peak.
            if peaks(1) <= x(k)
                indices(1) = [];
                peaks(1) = [];
                break;
            end
            
            % Advance.
            k = k - 1;
        end
        
        % If we have a peak at the end, check the wrapping portion just
        % before the start.
        k = 1;
        while chain_code_lengths(length(x)) - ...
                chain_code_lengths(indices(end)) + chain_code_lengths(k) < dist
            
            % Remove the peak.
            if peaks(end) < x(k)
                indices(end) = [];
                peaks(end) = [];
                break;
            end
            
            % Advance.
            k = k + 1;
        end
    end
end
end