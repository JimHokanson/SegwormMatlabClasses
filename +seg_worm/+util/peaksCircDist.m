function [peaks,indices] = peaksCircDist(x, dist,use_max,value_cutoff,chain_code_lengths)
%MAXPEAKSCIRCDIST Find the maximum peaks in a circular vector. The peaks
%are separated by, at least, the given distance.
%
%   [peaks,indices] = seg_worm.util.peaksCircDist(x, dist, *chain_code_lengths)
%
%       I think this algorithm grabs the best peaks, starting with the
%       first, then grabbing the 2nd (if we can), then the 3rd, etc ...
%
%       Importantly, a peak is not grabbed if it is not the max or min over
%       the window in which it is operating
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

if size(x,1) > 1
    x = x';
end


%Let's say we have the following distances, all relative to the first index
%-10 -3 -2 0 1 3 5 10 15 30 etc <- this is computed from x, see getLinearDistances
%          s
%          1 2 3 4 5  6  7 <- indices of the regular worm
%
%   Let's say our data length is 50 (length(x) => 50), then the following 
%   are indices on the front pad
%48  49 50 (corresponding to -10 -3 -2)
%
%   So, for example, if go +/- 6 on the first sampple, this would take us from
%   49 to 4 (distances -3 to 5). Using the indices output of 
%   getLinearDistances() we can grab 49,50,1,2,3,4 really easily. When
%   choosing index 1, we thus say we can no longer grab these other
%   indices.
%
%   NOTE: The distances are padded on both sides, not just the front


[distances,indices,~,x_locations] = seg_worm.util.getLinearDistances(chain_code_lengths);

x_cutoff_left  = x_locations - dist;
x_cutoff_right = x_locations + dist; 

n_distances = length(distances);
%round for both of these instead????
%It probably doesn't matter ...
start_I     = ceil(interp1(distances,1:n_distances,x_cutoff_left));
end_I       = floor(interp1(distances,1:n_distances,x_cutoff_right));

if use_max
   minMaxFH = @max;
else
   minMaxFH = @min;
end

n_sort = length(I);
n_x    = length(x);
is_peak_mask = false(1,n_x); 
taken_mask   = false(1,n_x); %Indicates that the index is close to or is 
%a peak and thus can not be used as a peak
for iElem = 1:n_sort
    cur_index = I(iElem);
    if ~taken_mask(cur_index)
       %NOTE: Even if this isn't the local max, it is greater
       %than anything that is by it, so it prevents anything
       %else from being used, so we might as well mark those indices
       %within it's distance as taken as well
       temp_indices             = indices(start_I(cur_index):end_I(cur_index));
       taken_mask(temp_indices) = true;
       
       %We need to ensure that within it's window of operatin that the
       %currently chosen value is a max or a min
       is_peak_mask(cur_index)  = minMaxFH(x(temp_indices)) == x(cur_index);
    end
end

indices = find(is_peak_mask);
peaks   = x(indices);

%This is only for max
%[peaks2,indices2] = helper__oldCodeMax(x, dist, chain_code_lengths,2*dist+1);

end

function [peaks,indices] = helper__oldCodeMax(x, dist, chain_code_lengths,winSize) %#ok<DEFNU>
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