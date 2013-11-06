function eventStats = removePartialEvents(eventStats, totalFrames)
%removePartialEvents  Remove partial events at the start and end of the data.
%
%   eventStats = removePartialEvents(eventStats, totalFrames)
%
%   Removes events that start right as the video starts or end right as the
%   video ends. Since the video does not extend further in either direction
%   we can't know for sure that the video has ended.
%
%
%   Called by:
%   seg_worm.w.stats.worm2histogram
%
%
%   seg_worm.feature.removePartialEvents
%
%   Inputs:
%       eventStats  - the event statistics (see events2stats)
%       totalFrames - the total number of frames in the video
%
%   Output:
%       eventStats - the event statistics with partial
%                    events, at the start and end of the data, removed
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

% Remove partially recorded events.
if ~isempty(eventStats)
    if eventStats(1).start == 0
        eventStats(1) = [];
    end
end
if ~isempty(eventStats)
    if eventStats(end).end == totalFrames
        eventStats(end) = [];
    end
end
end

