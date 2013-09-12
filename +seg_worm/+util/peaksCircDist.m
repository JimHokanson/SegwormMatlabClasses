function [peaks,indices] = peaksCircDist(x, dist,use_max,value_cutoff,chain_code_lengths)
%peaksCircDist  Find the maximum peaks in a circular vector. The peaks
%are separated by, at least, the given distance.
%
%   [peaks,indices] = seg_worm.util.peaksCircDist(x, dist, *chain_code_lengths)
%
%   This algorithm grabs the best peaks, starting with the first, then
%   grabbing the 2nd (if we can), then the 3rd, etc ...
%
%   Importantly, a peak is not grabbed if it is not the max or min over the
%   window in which it is operating
%
%   Inputs:
%       x                - the vector of values, I think these are
%                       generally angles
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

%The general algorithm is to start at the largest (or smallest if using
%min) and to progressively work down the scale

if ~exist('chain_code_lengths','var') || isempty(chain_code_lengths)
    % Use the array indices for length.
    chain_code_lengths = 1:length(x);
end

if size(x,1) > 1
    x  = x';
    is_transposed = true;
else
    is_transposed = false;
end

%xt - "x for testing" in some places in the code below
%it will be quicker (and/or easier) to assume that we want the largest
%value. By negating the data we can look for maxima (which will tell us
%where the minima are)
if ~use_max
    xt = -1*x;
else
    xt = x;
end

[start_I,end_I,indices] = helper__getExtentsAndIndices(chain_code_lengths,dist);

could_be_a_peak = helper__initCouldBeAPeak(min(end_I - start_I),xt);

if use_max
    I1 = find(x > value_cutoff & could_be_a_peak);
    [~,I2] = sort(-1*x(I1));
    I = I1(I2);
else
    I1     = find(x < value_cutoff & could_be_a_peak);
    [~,I2] = sort(x(I1));
    I = I1(I2);
end

is_peak_mask   = false(size(x));
%a peak and thus can not be used as a peak
n_sort = length(I);
for iElem = 1:n_sort
    cur_index = I(iElem);
    if could_be_a_peak(cur_index)
        %NOTE: Even if a point isn't the local max, it is greater
        %than anything that is by it that is currently not taken
        %(because of sorting), so it prevents these points
        %from undergoing the expensive search of determining
        %whether they are the min or max within their
        %else from being used, so we might as well mark those indices
        %within it's distance as taken as well
        temp_indices = indices(start_I(cur_index):end_I(cur_index));
        could_be_a_peak(temp_indices) = false;
        
        is_peak_mask(cur_index) = max(xt(temp_indices)) == xt(cur_index);
    end
end

indices = find(is_peak_mask);
peaks   = x(indices);

if is_transposed
    indices = indices';
    peaks = peaks';
end

%OLD DEBUGGING CODE
% [peaks2,indices2] = helper__oldCodeMax(x, dist, chain_code_lengths,2*dist+1);

end

function could_be_a_peak = helper__initCouldBeAPeak(min_dist,xt)
%We'll start with the extreme-most value and progress
%In the loop this allows us to set the could_be_a_peak to any
%of the neighbors as false because the value came first and is an extreme

%NOTE: We could make a loop of this ...
%We'll assume that we at least want to test the neighbor
%Corresponds to a min_dist of 3 (sort of, could be skewed)
%
%                   > right                   > left
%
could_be_a_peak = xt > [xt(2:end) xt(1)] & xt > [xt(end) xt(1:end-1)];

%This seems to be beneficial
if min_dist > 5
   %                                     > 2 to the right            > 2 to the left
   could_be_a_peak = could_be_a_peak & xt > [xt(3:end) xt(1:2)] & xt > [xt(end-1:end) xt(1:end-2)]; 
end


% %Doesn't pay off after the first few tests ...
% next_min    = 5;
% next_offset = 2;
% while next_min <= min_dist
%    could_be_a_peak = could_be_a_peak & xt > [xt(next_offset+1:end) xt(1:next_offset)] & xt > [xt(end-next_offset+1:end) xt(1:end-next_offset)]; 
%    next_offset = next_offset + 1;
%    next_min    = next_min + 2;
% end
%

% if min_dist > 7
%    %                                     > 3 to the right            > 3 to the left
%    could_be_a_peak = could_be_a_peak & xt > [xt(4:end) xt(1:3)] & xt > [xt(end-2:end) xt(1:end-3)]; 
% end
end
function [start_I,end_I,indices] = helper__getExtentsAndIndices(chain_code_lengths,dist)

persistent cc_input dist_input p_start p_end p_indices

if isequal(chain_code_lengths,cc_input) && isequal(dist,dist_input)
    indices = p_indices;
    start_I = p_start;
    end_I   = p_end;
else
    
    [distances,indices,~,x_locations] = seg_worm.util.getLinearDistances(chain_code_lengths);
    
    x_cutoff_left  = x_locations - dist;
    x_cutoff_right = x_locations + dist;
    
    n_distances = length(distances);
    %round for both of these instead????
    %It probably doesn't matter ...
    
    %For every point, we find the indices over which it must be the maximum (or
    %minimum). This is based on distances, not indices
    F = griddedInterpolant(distances,1:n_distances);
    
    start_I = ceil(F(x_cutoff_left));
    end_I   = floor(F(x_cutoff_right));
    
    cc_input   = chain_code_lengths;
    dist_input = dist;
    p_start    = start_I;
    p_end      = end_I;
    p_indices  = indices;
    
end

end

function [peaks,indices] = helper__oldCodeMax(x, dist, chain_code_lengths,winSize)
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