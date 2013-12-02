function frames = findEvent(data, minThr, maxThr, varargin)
%FINDEVENT Find an event within an array of data.
%
%
%   Note: if the first/last event are solely preceded/followed by NaN
%   frames, these frames are swallowed into the respective event.
%
%   frames = seg_worm.feature.event.findEvent(data, minThr, maxThr, varargin);
%
%
%
%   ******************************************************
%       YIKES! This needs to be cleaned up ...
%   ******************************************************
%
%
%
%
%
%   Rows are functional groups of inputs ...
%
%   FRAMES = FINDEVENT(DATA, MINTHR, MAXTHR, 
%
%                      ISATTHR,
%
%                      MINFRAMESTHR, MAXFRAMESTHR, ISATFRAMESTHR,
%
%                      MINSUMTHR, MAXSUMTHR, ISATSUMTHR, 
%
%                      SUMDATA,
%
%                      MININTERFRAMESTHR, MAXINTERFRAMESTHR, ISATINTERFRAMESTHR,
%
%                      MININTERSUMTHR, MAXINTERSUMTHR, ISATINTERSUMTHR)
%
%
%   FRAMES = FINDEVENT(
%       DATA                1
%       MINTHR              2
%       MAXTHR,             3
%       ISATTHR,            4
%       MINFRAMESTHR,       5
%       MAXFRAMESTHR,       6
%       ISATFRAMESTHR,      7
%       MINSUMTHR,          8
%       MAXSUMTHR,          9
%       ISATSUMTHR,         10
%       SUMDATA,            11
%       MININTERFRAMESTHR,  12
%       MAXINTERFRAMESTHR,  13
%       ISATINTERFRAMESTHR, 14
%       MININTERSUMTHR,     15
%       MAXINTERSUMTHR,     16
%       ISATINTERSUMTHR)    17
%
%
%
%
%
%   Inputs:
%       data              - the event data
%       minThr            - the minimum threshold for the data;
%                           if empty, their is no minimum threshold
%       maxThr            - the maximum threshold for the data
%                           if empty, their is no maximum threshold
%
%--------------------------------------------------------------------------
%
%       isAtThr           - (false) is the data threshold inclusive?
%
%--------------------------------------------------------------------------
%
%       minFramesThr      - ([], no minimum) the minimum threshold for the event frames
%       maxFramesThr      - ([], no maximum) the maximum threshold for the event frames
%       isAtFramesThr     - (false) is the frames threshold inclusive?
%
%--------------------------------------------------------------------------
%
%       minSumThr         - ([], no minimum) the minimum threshold for the data sum
%       maxSumThr         - ([], no maximum) the maximum threshold for the data sum
%       isAtSumThr        - (false) is the sum threshold inclusive?
%
%--------------------------------------------------------------------------
%
%       sumData           - the event data for sum thresholding;
%                           if empty, the event data is used
%
%--------------------------------------------------------------------------
%
%       minInterFramesThr - the minimum threshold for the number of frames
%                           separating events, events separated by less
%                           frames are unified;
%                           if empty, their is no minimum threshold
%       maxInterFramesThr - the maximum threshold for the number of frames
%                           separating events, events separated by more
%                           frames are unified;
%                           if empty, their is no maximum threshold
%       isAtInterFramesThr - (false) is the frame separation threshold inclusive?
%
%--------------------------------------------------------------------------
%
%
%       minInterSumThr    - the minimum threshold for the data sum
%                           separating events, events separated by less
%                           data are unified;
%                           if empty, their is no minimum threshold
%       maxInterSumThr    - the maximum threshold for the data sum
%                           separating events, events separated by more
%                           data are unified;
%                           if empty, their is no minimum threshold
%       isAtInterSumThr   - (false) is the data sum separation threshold inclusive?
%
%--------------------------------------------------------------------------
%   Output:
%       frames - the frames at which the event took place;
%                a structure array with fields:
%
%                start = the start frame
%                end   = the end frame
%
%
%   Called by:
%   seg_worm.feature_helpers.locomotion.getWormMotionCodes
%
%   See also:
%   EVENTS2STATS, EVENTS2ARRAY
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software;
% you must reproduce all copyright notices and other proprietary
% notices on any copies of the Software.

%TODO: Who calls this function?
%seg_worm.feature_helpers.locomotion.getWormMotionCodes
%
%call form:
%forwardFrames ... (midbody_speed, minForwardSpeed, [], 
    %true, 
    %wormEventFramesThr, 
    %[], 
    %false, ...
    %minForwardDistance, 
    %[], 
    %true, 
    %distance, 
    %wormEventMinInterFramesThr);


%This is for verification that we're doing things right, will be removed
%...
frames2 = h__findEvent_old(data, minThr, maxThr, varargin{:});

%==========================================================================
% Is the threshold inclusive?
isAtThr = false;
if ~isempty(varargin)
    isAtThr = varargin{1};
end

% Determine the minimum event frames threshold.
minFramesThr = [];
if length(varargin) > 1
    minFramesThr = varargin{2}(:);
end

% Determine the maximum event frames threshold.
maxFramesThr = [];
if length(varargin) > 2
    maxFramesThr = varargin{3}(:);
end

% Is the event frames threshold inclusive?
isAtFramesThr = false;
if length(varargin) > 3
    isAtFramesThr = varargin{4};
end

% Determine the minimum data sum threshold.
minSumThr = [];
if length(varargin) > 4
    minSumThr = varargin{5}(:);
end

% Determine the maximum data sum threshold.
maxSumThr = [];
if length(varargin) > 5
    maxSumThr = varargin{6}(:);
end

% Is the data sum threshold inclusive?
isAtSumThr = false;
if length(varargin) > 6
    isAtSumThr = varargin{7}(:);
end

% Determine the event data for sum thresholding.
sumData = [];
if length(varargin) > 7
    sumData = varargin{8}(:);
end
if isempty(sumData)
    sumData = data;
end

% Determine the minimum frames separation threshold.
minInterFramesThr = [];
if length(varargin) > 8
    minInterFramesThr = varargin{9}(:);
end

% Determine the maximum frames separation threshold.
maxInterFramesThr = [];
if length(varargin) > 9
    maxInterFramesThr = varargin{10}(:);
end

% Is the frames separation threshold inclusive?
isAtInterFramesThr = false;
if length(varargin) > 10
    isAtInterFramesThr = varargin{11};
end

% Determine the minimum data sum separation threshold.
minInterSumThr = [];
if length(varargin) > 11
    minInterSumThr = varargin{12}(:);
end

% Determine the maximum data sum separation threshold.
maxInterSumThr = [];
if length(varargin) > 12
    maxInterSumThr = varargin{13}(:);
end

% Is the data sum separation threshold inclusive?
isAtInterSumThr = false;
if length(varargin) > 13
    isAtInterSumThr = varargin{14};
end
%==========================================================================


% Fix the data.
%--------------------------------------------------------------------------
data    = data(:);
sumData = sumData(:);
minThr  = minThr(:);
maxThr  = maxThr(:);

% Mark the events.
%--------------------------------------------------------------------------
isEvent = true(length(data),1);
if ~isempty(minThr)
    if isAtThr
        isEvent = data >= minThr;
    else
        isEvent = data > minThr;
    end
end

if ~isempty(maxThr)
    if isAtThr
        isEvent = isEvent & data <= maxThr;
    else
        isEvent = isEvent & data < maxThr;
    end
end

% Find the events.
%----------------------------------------------------
dEvent      = diff(isEvent);
startFrames = find(dEvent == 1);
endFrames   = find(dEvent == -1) - 1;

if isEvent(1)
    startFrames(2:(end + 1)) = startFrames;
    startFrames(1) = 1;
end
if isEvent(end)
    endFrames(end + 1) = length(isEvent);
end


%Possible short circuit ...
%--------------------------------------------------------------------------
if isempty(startFrames)
    frames = [];
    return
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Do we have any events?

% Include NaNs at the start and end.
%----------------------------------------------------------------------
if startFrames(1) > 0 && all(isnan(data(1:startFrames(1))))
    startFrames(1) = 1;
end

if endFrames(end) < length(isEvent) - 1 && all(isnan(data((endFrames(end) + 2):end)))
    endFrames(end) = length(isEvent);
end

%--------------------------------------------------------------------------
[startFrames,endFrames] = h__unifyEvents(startFrames,endFrames,minInterFramesThr,maxInterFramesThr,isAtInterFramesThr);


%--------------------------------------------------------------------------
%Is this really the same thing twice with different values ????
%I'm  99% sure this isn't done right
if ~isempty(minInterSumThr) || ~isempty(maxInterSumThr)
    error('I don''t think this was coded right to start ..., check code')
end
%NOTE: this should use data, but it doesn't
%Perhaps there exists a correct version?
[startFrames,endFrames] = h__unifyEvents(startFrames,endFrames,minInterSumThr,maxInterSumThr,isAtInterSumThr);


%--------------------------------------------------------------------------
[startFrames,endFrames] = h__removeTooSmallOrLargeEvents(startFrames,endFrames,minFramesThr,maxFramesThr,isAtFramesThr);




% Check the event sums.
%--------------------------------------------------------------------------
if ~(isempty(minSumThr) && isempty(maxSumThr))
    
    % Compute the event sums.
    eventSums = nan(length(startFrames), 1);
    for i = 1:length(eventSums)
        eventSums(i) = nansum(sumData((startFrames(i)):(endFrames(i))));
    end
    
    % Compute the event sum thresholds.
    if length(minSumThr) > 1
        newMinSumThr = nan(size(eventSums));
        for i = 1:length(newMinSumThr)
            newMinSumThr(i) = nanmean(minSumThr((startFrames(i)):(endFrames(i))));
        end
        minSumThr = newMinSumThr;
    end
    
    if length(maxSumThr) > 1
        newMaxSumThr = nan(size(eventSums));
        for i = 1:length(newMaxSumThr)
            newMaxSumThr(i) = nanmean(maxSumThr((startFrames(i)):(endFrames(i))));
        end
        maxSumThr = newMaxSumThr;
    end
    
    % Remove small events.
    removeEvents = false(size(eventSums));
    if ~isempty(minSumThr)
        if isAtSumThr
            removeEvents = eventSums <= minSumThr;
        else
            removeEvents = eventSums < minSumThr;
        end
    end
    
    % Remove large events.
    if ~isempty(maxSumThr)
        if isAtSumThr
            removeEvents =  removeEvents | eventSums >= maxSumThr;
        else
            removeEvents =  removeEvents | eventSums > maxSumThr;
        end
    end
    
    % Remove the events.
    startFrames(removeEvents) = [];
    endFrames(removeEvents)   = [];
end  %~(isempty(minSumThr) && isempty(maxSumThr))

if isempty(startFrames)
    frames = [];
else
    frames = struct('start',num2cell(startFrames),'end',num2cell(endFrames));
end

end

%==========================================================================







%==========================================================================


function [startFrames,endFrames] = h__unifyEvents(startFrames,endFrames,minInterFramesThr,maxInterFramesThr,isAtInterFramesThr)

% Check the number of frames separating events.
% Note: unify events before checking their frames and sums.

% Unify small time gaps.
if ~isempty(minInterFramesThr)
    if isAtInterFramesThr
        [startFrames,endFrames] = h__removeGaps(startFrames,endFrames,minInterFramesThr,@le); %  <=
    else % the threshold is exclusive
        [startFrames,endFrames] = h__removeGaps(startFrames,endFrames,minInterFramesThr,@lt); %  <
    end
end

% Unify large time gaps.
if ~isempty(maxInterFramesThr)
    if isAtInterFramesThr
        [startFrames,endFrames] = h__removeGaps(startFrames,endFrames,maxInterFramesThr,@ge); %  >=
    else
        % the threshold is exclusive
        [startFrames,endFrames] = h__removeGaps(startFrames,endFrames,maxInterFramesThr,@gt); %  >=
    end
end

end


function [startFrames,endFrames] = h__removeGaps(startFrames,endFrames,right_comparison_value,fh)

% Find small gaps.
i = 1;
while i < length(startFrames)
    
    % Swallow the gaps.
    %NOTE: fh is <, <=, >, >=
    while i < length(startFrames) && fh(startFrames(i + 1) - endFrames(i) - 1,right_comparison_value)
        endFrames(i)       = endFrames(i + 1);
        startFrames(i + 1) = [];
        endFrames(i + 1)   = [];
    end
    
    % Advance.
    i = i + 1;
end
end

function [startFrames,endFrames] = h__removeTooSmallOrLargeEvents(startFrames,endFrames,minFramesThr,maxFramesThr,isAtFramesThr)

% Check the event frames.
%--------------------------------------------------------------------------
if ~(isempty(minFramesThr) && isempty(maxFramesThr))
    
    % Compute the event frames.
    n_frames_per_event = endFrames - startFrames + 1;
    
    % Remove small events.
    removeEvents = false(size(n_frames_per_event));
    if ~isempty(minFramesThr)
        if isAtFramesThr
            removeEvents = n_frames_per_event <= minFramesThr;
        else
            removeEvents = n_frames_per_event < minFramesThr;
        end
    end
    
    % Remove large events.
    if ~isempty(maxFramesThr)
        if isAtFramesThr
            removeEvents = removeEvents | n_frames_per_event >= maxFramesThr;
        else
            removeEvents = removeEvents | n_frames_per_event > maxFramesThr;
        end
    end
    
    % Remove the events.
    startFrames(removeEvents) = [];
    endFrames(removeEvents)   = [];
end

end








































































































function frames = h__findEvent_old(data, minThr, maxThr, varargin)

%==========================================================================
% Is the threshold inclusive?
isAtThr = false;
if ~isempty(varargin)
    isAtThr = varargin{1};
end

% Determine the minimum event frames threshold.
minFramesThr = [];
if length(varargin) > 1
    minFramesThr = varargin{2}(:);
end

% Determine the maximum event frames threshold.
maxFramesThr = [];
if length(varargin) > 2
    maxFramesThr = varargin{3}(:);
end

% Is the event frames threshold inclusive?
isAtFramesThr = false;
if length(varargin) > 3
    isAtFramesThr = varargin{4};
end

% Determine the minimum data sum threshold.
minSumThr = [];
if length(varargin) > 4
    minSumThr = varargin{5}(:);
end

% Determine the maximum data sum threshold.
maxSumThr = [];
if length(varargin) > 5
    maxSumThr = varargin{6}(:);
end

% Is the data sum threshold inclusive?
isAtSumThr = false;
if length(varargin) > 6
    isAtSumThr = varargin{7}(:);
end

% Determine the event data for sum thresholding.
sumData = [];
if length(varargin) > 7
    sumData = varargin{8}(:);
end
if isempty(sumData)
    sumData = data;
end

% Determine the minimum frames separation threshold.
minInterFramesThr = [];
if length(varargin) > 8
    minInterFramesThr = varargin{9}(:);
end

% Determine the maximum frames separation threshold.
maxInterFramesThr = [];
if length(varargin) > 9
    maxInterFramesThr = varargin{10}(:);
end

% Is the frames separation threshold inclusive?
isAtInterFramesThr = false;
if length(varargin) > 10
    isAtInterFramesThr = varargin{11};
end

% Determine the minimum data sum separation threshold.
minInterSumThr = [];
if length(varargin) > 11
    minInterSumThr = varargin{12}(:);
end

% Determine the maximum data sum separation threshold.
maxInterSumThr = [];
if length(varargin) > 12
    maxInterSumThr = varargin{13}(:);
end

% Is the data sum separation threshold inclusive?
isAtInterSumThr = false;
if length(varargin) > 13
    isAtInterSumThr = varargin{14};
end
%==========================================================================


% Fix the data.
%--------------------------------------------------------------------------
data    = data(:);
sumData = sumData(:);
minThr  = minThr(:);
maxThr  = maxThr(:);

% Mark the events.
%--------------------------------------------------------------------------
isEvent = true(length(data),1);
if ~isempty(minThr)
    if isAtThr
        isEvent = data >= minThr;
    else
        isEvent = data > minThr;
    end
end

if ~isempty(maxThr)
    if isAtThr
        isEvent = isEvent & data <= maxThr;
    else
        isEvent = isEvent & data < maxThr;
    end
end

% Find the events.
dEvent      = diff(isEvent);
startFrames = find(dEvent == 1);
endFrames   = find(dEvent == -1) - 1;

% All event must have a start and end.
% if isempty(startFrames)
%     endFrames = [];
% elseif isempty(endFrames)
%     startFrames= [];
% else
%     if startFrames(1) > endFrames(1)
%         endFrames(1) = [];
%     end
%     if startFrames(end) > endFrames(end)
%         startFrames(end) = [];
%     end
% end

if isEvent(1)
    startFrames(2:(end + 1)) = startFrames;
    startFrames(1) = 0;
end
if isEvent(end)
    endFrames(end + 1) = length(isEvent) - 1;
end


%Possible short circuit ...
%--------------------------------------------------------------------------
if isempty(startFrames)
    frames = [];
    return
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Do we have any events?

% Include NaNs at the start and end.
%----------------------------------------------------------------------
if startFrames(1) > 0 && all(isnan(data(1:startFrames(1))))
    startFrames(1) = 0;
end

if endFrames(end) < length(isEvent) - 1 && all(isnan(data((endFrames(end) + 2):end)))
    endFrames(end) = length(isEvent) - 1;
end

% Check the number of frames separating events.
% Note: unify events before checking their frames and sums.
if ~(isempty(minInterFramesThr) && isempty(maxInterFramesThr))
    
    % Unify small time gaps.
    if ~isempty(minInterFramesThr)
        if isAtInterFramesThr
            
            % Find small gaps.
            i = 1;
            while i < length(startFrames)
                
                % Swallow the gaps.
                while i < length(startFrames) && startFrames(i + 1) - endFrames(i) - 1 <= minInterFramesThr
                    endFrames(i) = endFrames(i + 1);
                    startFrames(i + 1) = [];
                    endFrames(i + 1) = [];
                end
                
                % Advance.
                i = i + 1;
            end
            
        else % the threshold is exclusive
            
            % Find small gaps.
            i = 1;
            while i < length(startFrames)
                
                % Swallow the gaps.
                while i < length(startFrames) && startFrames(i + 1) - endFrames(i) - 1 < minInterFramesThr
                    endFrames(i)       = endFrames(i + 1);
                    startFrames(i + 1) = [];
                    endFrames(i + 1)   = [];
                end
                
                % Advance.
                i = i + 1;
            end
        end
    end
    
    % Unify large time gaps.
    if ~isempty(maxInterFramesThr)
        if isAtInterFramesThr
            
            % Find large gaps.
            i = 1;
            while i < length(startFrames)
                
                % Swallow the gaps.
                while i < length(startFrames) && startFrames(i + 1) - endFrames(i) - 1 >= maxInterFramesThr
                    endFrames(i)       = endFrames(i + 1);
                    startFrames(i + 1) = [];
                    endFrames(i + 1)   = [];
                end
                
                % Advance.
                i = i + 1;
            end
            
        else % the threshold is exclusive
            
            % Find large gaps.
            i = 1;
            while i < length(startFrames)
                
                % Swallow the gaps.
                while i < length(startFrames) && startFrames(i + 1) - endFrames(i) - 1 > maxInterFramesThr
                    endFrames(i)       = endFrames(i + 1);
                    startFrames(i + 1) = [];
                    endFrames(i + 1)   = [];
                end
                
                % Advance.
                i = i + 1;
            end
        end
    end
end

% Check the data sum separating events.
% Note: unify events before checking their frames and sums.
if ~(isempty(minInterSumThr) && isempty(maxInterSumThr))
    
    % Unify small time gaps.
    if ~isempty(minInterSumThr)
        if isAtInterSumThr
            
            % Find small gaps.
            i = 1;
            while i < length(startFrames)
                
                % Swallow the gaps.
                while i < length(startFrames) && ...
                        startFrames(i + 1) - endFrames(i) - 1 <= ...
                        minInterSumThr
                    endFrames(i) = endFrames(i + 1);
                    startFrames(i + 1) = [];
                    endFrames(i + 1) = [];
                end
                
                % Advance.
                i = i + 1;
            end
            
        else % the threshold is exclusive
            
            % Find small gaps.
            i = 1;
            while i < length(startFrames)
                
                % Swallow the gaps.
                while i < length(startFrames) && startFrames(i + 1) - endFrames(i) - 1 < minInterSumThr
                    endFrames(i)       = endFrames(i + 1);
                    startFrames(i + 1) = [];
                    endFrames(i + 1)   = [];
                end
                
                % Advance.
                i = i + 1;
            end
        end
    end
    
    % Unify large time gaps.
    if ~isempty(maxInterSumThr)
        if isAtInterSumThr
            
            % Find large gaps.
            i = 1;
            while i < length(startFrames)
                
                % Swallow the gaps.
                while i < length(startFrames) && startFrames(i + 1) - endFrames(i) - 1 >= maxInterSumThr
                    endFrames(i) = endFrames(i + 1);
                    startFrames(i + 1) = [];
                    endFrames(i + 1) = [];
                end
                
                % Advance.
                i = i + 1;
            end
            
        else % the threshold is exclusive
            
            % Find large gaps.
            i = 1;
            while i < length(startFrames)
                
                % Swallow the gaps.
                while i < length(startFrames) && startFrames(i + 1) - endFrames(i) - 1 > maxInterSumThr
                    endFrames(i) = endFrames(i + 1);
                    startFrames(i + 1) = [];
                    endFrames(i + 1) = [];
                end
                
                % Advance.
                i = i + 1;
            end
        end
    end
end

% Check the event frames.
%--------------------------------------------------------------------------
if ~(isempty(minFramesThr) && isempty(maxFramesThr))
    
    % Compute the event frames.
    eventNumFrames = endFrames - startFrames + 1;
    
    % Remove small events.
    removeEvents = false(size(eventNumFrames));
    if ~isempty(minFramesThr)
        if isAtFramesThr
            removeEvents = eventNumFrames <= minFramesThr;
        else
            removeEvents = eventNumFrames < minFramesThr;
        end
    end
    
    % Remove large events.
    if ~isempty(maxFramesThr)
        if isAtFramesThr
            removeEvents = removeEvents || eventNumFrames >= maxFramesThr;
        else
            removeEvents = removeEvents || eventNumFrames > maxFramesThr;
        end
    end
    
    % Remove the events.
    startFrames(removeEvents) = [];
    endFrames(removeEvents) = [];
end

% Check the event sums.
%--------------------------------------------------------------------------
if ~(isempty(minSumThr) && isempty(maxSumThr))
    
    % Compute the event sums.
    eventSums = nan(length(startFrames), 1);
    for i = 1:length(eventSums)
        eventSums(i) = nansum(sumData((startFrames(i) + 1):(endFrames(i) + 1)));
    end
    
    % Compute the event sum thresholds.
    if length(minSumThr) > 1
        newMinSumThr = nan(size(eventSums));
        for i = 1:length(newMinSumThr)
            newMinSumThr(i) = nanmean(minSumThr((startFrames(i) + 1):(endFrames(i) + 1)));
        end
        minSumThr = newMinSumThr;
    end
    if length(maxSumThr) > 1
        newMaxSumThr = nan(size(eventSums));
        for i = 1:length(newMaxSumThr)
            newMaxSumThr(i) = nanmean(maxSumThr((startFrames(i) + 1):(endFrames(i) + 1)));
        end
        maxSumThr = newMaxSumThr;
    end
    
    % Remove small events.
    removeEvents = false(size(eventSums));
    if ~isempty(minSumThr)
        if isAtSumThr
            removeEvents = eventSums <= minSumThr;
        else
            removeEvents = eventSums < minSumThr;
        end
    end
    
    % Remove large events.
    if ~isempty(maxSumThr)
        if isAtSumThr
            removeEvents =  removeEvents || eventSums >= maxSumThr;
        else
            removeEvents =  removeEvents || eventSums > maxSumThr;
        end
    end
    
    % Remove the events.
    startFrames(removeEvents) = [];
    endFrames(removeEvents)   = [];
end  %~(isempty(minSumThr) && isempty(maxSumThr))

% Record the event frames.
frames = struct( ...
    'start', [], ...
    'end',  []);
for i = 1:length(startFrames)
    frames(i).start = startFrames(i);
    frames(i).end = endFrames(i);
end

% Did we find any event frames?
if isempty(frames(1).start)
    frames = [];
end


end
