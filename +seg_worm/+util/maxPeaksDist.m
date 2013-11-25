function [peaks, indices] = maxPeaksDist(x, dist,use_max,value_cutoff,chain_code_lengths)
%MAXPEAKSDIST Find the maximum peaks in a vector. The peaks
%are separated by, at least, the given distance.
%
%   seg_worm.util.maxPeaksDist
%   
%   [PEAKS INDICES] = seg_worm.util.maxPeaksDist(x, dist,use_max,value_cutoff,*chain_code_lengths)
%
%   Inputs:
%       x                - the vector of values
%       dist             - the minimum distance between peaks
%       chainCodeLengths - the chain-code length at each index;
%                          if empty, the array indices are used instead
%
%   Outputs:
%       peaks   - the maximum peaks
%       indices - the indices for the peaks
%
%   See also MINPEAKSDIST, COMPUTECHAINCODELENGTHS
%
%   ****************
%   Used in seg_worm.feature_helpers.posture.getAmplitudeAndWavelength
%
%   NOTE: Outputs ARE NOT SORTED 
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

%JAH TODO: Merge with circular version ...
%Create shared functions ...


%   ????? - are the outputs sorted ?????? i.e. peaks(1) > peaks(2) >
%   peaks(3) ????
%
%   I don't think so ....
%   
%{
  x = zeros(1,100);
  x(1) = 10;
  x(10) = 20;
  x(30) = 30;
  x(50) = 15;
  dist = 5;
  [peaks, indices] = seg_worm.util.maxPeaksDist(x, dist);
%}


%NOTE: This is very similar to:
%seg_worm.util.peaksCircDist


%JAH: I think the original implementation of this code had an error in it
%
%JAH: Speed note, I think it would be faster when indices are used
%to run a filter on neighbors +/- 1, I tried this in the peaksCircDist
%and it helped a little bit but I removed it because it complicated the
%code

% Are there chain-code lengths?
if ~exist('chain_code_lengths','var') || isempty(chain_code_lengths)
    chain_code_lengths = 1:length(x);
else
    %start and ends need to be adjusted ...
    error('Code as implemented is not correct for this ...')
end

% Is the vector larger than the search window?
winSize = 2 * dist + 1;
if chain_code_lengths(end) < winSize
    [peaks, indices] = max(x);
    return;
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

if use_max
    could_be_a_peak = x > value_cutoff;
    I1 = find(could_be_a_peak);
    [~,I2] = sort(-1*x(I1));
    I = I1(I2);
else
    could_be_a_peak = x < value_cutoff;
    I1     = find(could_be_a_peak);
    [~,I2] = sort(x(I1));
    I = I1(I2);
end

n_points = length(x);

%This code would need to be fixed if real distances
%are input ...
too_close = dist - 1;

temp_I  = 1:n_points;
start_I = temp_I - too_close; %Note, separated by dist is ok
%This sets the off limits area, so we go in by 1
end_I   = temp_I + too_close;


start_I(start_I < 1) = 1;
end_I(end_I > n_points) = n_points;

is_peak_mask   = false(1,n_points);
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
        temp_indices = start_I(cur_index):end_I(cur_index);
        could_be_a_peak(temp_indices) = false;
        
        is_peak_mask(cur_index) = max(xt(temp_indices)) == xt(cur_index);
    end
end

indices = find(is_peak_mask);
peaks   = x(indices);

% % % [peaks2, indices2] = h_old(x,dist,chain_code_lengths);
% % % 
% % % mask = peaks2 > 0.5*max(peaks2);
% % % indices2 = indices2(mask);
% % % peaks2 = peaks2(mask);
% % % 
% % % if ~isequal(peaks(:),peaks2(:)) || ~isequal(indices(:),indices2(:))
% % %     error('wtf')
% % % end

% keyboard

end


function [peaks, indices] = h_old(x,dist,chainCodeLengths)

%NOTE: I think there is a off by 1 bug here ...


% Initialize the peaks and indices.
winSize = 2 * dist + 1;
wins    = ceil(length(x) / winSize);
peaks   = zeros(wins, 1); % pre-allocate memory
indices = zeros(wins, 1); % pre-allocate memory

% Search for peaks.
im = 0; % the last maxima index
ie = 0; % the end index for the last maxima's search window
ip = 1; % the current, potential, max peak index
p = x(ip); % the current, potential, max peak value
i = 2; % the vector index
j = 1; % the recorded, maximal peaks index
while i <= length(x)
    
    % Found a potential peak.
    if (isnan(p) && ~isnan(x(i))) || x(i) > p
        ip = i;
        p = x(i);
    end
    
    % Test the potential peak.
    if ~isnan(p) && (chainCodeLengths(i) - chainCodeLengths(ip) >= dist ...
            || i == length(x))
        
        % Check the untested values next to the previous maxima.
        if im > 0 && chainCodeLengths(ip) - chainCodeLengths(im) <= 2 * dist
            
            % Check the untested values next to the previous maxima. 
            isMax = true;
            k = ie;
            while isMax && k > 0 && chainCodeLengths(ip) - chainCodeLengths(k) < dist
                
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
end