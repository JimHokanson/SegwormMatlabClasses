function [eventStats, summaryStats] = events2stats(event_ss, fps, data, name, interName)
%EVENTS2STATS Compute the event statistics.
%
%   [eventStats,summaryStats] = seg_worm.events.events2stats(event_ss, fps, data, name, interName)
%
%
%
%   Inputs:
%       event_ss  - the event frames (see findEvent)
%                    .start
%                    .end
%   
%       fps       - the video's frames/second
%       data      - the data values
%       name      - the struct field name for the event's data sum;
%                   if empty, this value is neither computed nor included
%       interName - the struct field name for the data sum till the next event;
%                   if empty, this value is neither computed nor included
%
%   Outputs:
%       eventStats   - the event statistics; a structure array with fields:
%
%                      start       = the start frame
%                      end         = the end frame
%                      time        = the event time
%                      <name>      = the sum of the event data
%                      interTime   = the time till the next event
%                      inter<name> = the sum of the data till the next event
%
%       summaryStats - the summary statistics for the events;
%                      a structure with fields:
%
%                      frequency = the event frequency (excluding partial
%                                  events at the start and end)
%                      ratio     = a structure with fields:
%
%                                  time   = the ratio of time
%                                           (event time / total time)
%                                  <name> = the ratio of data
%                                           (event data / total data)
%
%   See also:
%   FINDEVENT
%
%   EXAMPLE CALLS:
%   %coiled statistics
%   [coilEventStats, coiledStats] = events2stats(coilFrames, fps, distance, [], 'interDistance');
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software;
% you must reproduce all copyright notices and other proprietary
% notices on any copies of the Software.

n_events     = length(event_ss);
n_vid_frames = length(data);

% Are there any event frames?
if n_events == 0
    
    eventStats = [];
    % There are no events.
    frequency = 0;
    if isName
        ratios = struct('time', 0, name, 0);
    else
        ratios = struct('time', 0);
    end
    
    % Create the summary statistics.
    summaryStats = struct('frequency', frequency,'ratio', ratios);
    return;
end

%--------------------------------------------------------------------------
% Fix the data.
data = data(:);

start_Is = [event_ss.start];
end_Is   = [event_ss.end];

elapsed_times = (end_Is - start_Is + 1)./fps;
inter_times   = [end_Is(2:end) - start_Is(1:end-1) - 1 NaN]./fps;

eventStats = struct( ...
    'start',        num2cell(start_Is), ...
    'end',          num2cell(end_Is), ...
    'time',         num2cell(elapsed_times), ...
    'interTime',    num2cell(inter_times));

%Data statistics if requested
%--------------------------------------------------------------------------
if ~isempty(name)
   for iEvent = 1:n_events
      eventStats(iEvent).(name) = nansum(data(start_Is(iEvent):end_Is(iEvent))); 
   end
end

if ~isempty(interName)
   for iEvent = 1:(n_events-1)
      start_frame = end_Is(iEvent)+1;
      end_frame   = start_Is(iEvent+1)-1;
      eventStats(iEvent).(interName) = nansum(data(start_frame:end_frame));  
   end
   eventStats(end).(interName) = NaN;
end


% Compute the number of events, excluding the partially recorded ones.

n_events_for_stats = n_events;

if n_events_for_stats > 1
    if start_Is(1) == 1
        n_events_for_stats = n_events_for_stats - 1;
    end
    if end_Is(end) == n_vid_frames
        n_events_for_stats = n_events_for_stats - 1;
    end
end

% Compute the event statistics.
totalTime = n_vid_frames / fps;
frequency = n_events_for_stats / totalTime;
timeRatio = nansum([eventStats.time]) / totalTime;

ratios.time = timeRatio;
if ~isempty(name)
   ratios.(name) = nansum([eventStats.(name)])/nansum(data); 
end

summaryStats = struct( ...
    'frequency', frequency, ...
    'ratio',     ratios);
end

