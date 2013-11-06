function [isUsed,wormX,tickX,nameX,labelX] = wormInfo2MonthOfYear(info, filt)
%wormInfo2MonthOfYear Convert worm information to the month the experiment
%was performed.
%
%   [isUsed,wormX,tickX,nameX,labelX] = wormInfo2MonthOfYear(info, filt)
%
%   Inputs:
%       info - the worm information
%       filt - the filtering criteria; a structure with any of the fields:
%              minFPS     = the minimum video frame rate (frames/seconds)
%              minTime    = the minimum video time (seconds)
%              maxTime    = the maximum video time (seconds)
%              minSegTime = the minimum time for segmented video (seconds)
%              minRatio   = the minimum ratio for segmented video frames
%              minDate    = the minimum date to use (DATENUM)
%              maxDate    = the maximum date to use (DATENUM)
%              years      = the years to use
%              months     = the months to use (1-12)
%              weeks      = the weeks to use (1-52)
%              days       = the days (of the week) to use (1-7)
%              hours      = the hours to use (1-24)
%              trackers   = the trackers to use (1-8)
%              ticks      = the number of ticks to use
%
%   Outputs:
%       isUsed - for each worm, are we using the data?
%       wormX  - the x-axis values
%       tickX  - the x-axis ticks
%       nameX  - the x-axis tick names
%       labelX - the x-axis tick labels
%
%   See also:
%   PLOTWORMSTATS, FILTERWORMINFO
%
%
% � Medical Research Council 2012
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

% Filter the worms.
[isUsed,data] = filterWormInfo(info, filt);

% Compute the x-axis values.
wormX(1,:) = [data(isUsed).month];

% Compute the x-axis ticks.
tickX = sort(unique(wormX));
if ~isempty(filt) && isstruct(filt)
    if isfield(filt, 'ticks')
        if ~isempty(filt.ticks) && filt.ticks < length(tickX)
            tickX = round(linspace(tickX(1), tickX(end), filt.ticks));
        end
    end
end

% Construct the x-axis tick names.
nameX(1,:) = arrayfun(@(x) [datestr(datenum(2000, x, 15), 'mmm') ' (N = ' ...
    num2str(sum(x == wormX)) ')'], tickX, 'UniformOutput', false);

% Construct the x-axis label.
labelX = [];
end