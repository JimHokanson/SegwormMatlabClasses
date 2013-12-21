function initObject(obj,data)
%
%   seg_worm.stats.hist.initObject()
%
%
%   This is essentially the constructor code. I moved it in here to avoid
%   the indenting.
%


%seg_worm.util.histogram
%seg_worm.util.nanHistogram

data(isinf(data)) = NaN;

n_samples = sum(~isnan(data));

obj.samples = n_samples;

if n_samples == 0
    return
end

[obj.bins,edges] = h_computeBinInfo(data,obj.resolution);

% Compute the histogram counts for all the data.
counts = histc(data, edges);
counts(end) = []; %Remove the extra bin at the end (for overflow)

obj.counts = counts;

obj.pdf = counts./sum(counts);

obj.mean = nanmean(data);
obj.std  = nanstd(data);

end

function [bins,edges] = h_computeBinInfo(data,resolution)
%
%
%   NOTE: This version may have an extra bin than the previous version but
%   this one is MUCH simpler and merging should be much simpler as edges
%   should always align ...
%
%   
%   min -65.4
%   max 20.01
%   resolution 1
%   Old:
%   edges -65.5 to 20.5
%   New:
%   edges -70 to 21
%
%   


MAX_N = 1e6;

%JAH: I find this code to be very confusing ...


%Compute the data range & padding
%------------------------------------------------
min_data = min(data);
max_data = max(data);

min_edge = floor(min_data/resolution)*resolution;
max_edge = ceil(max_data/resolution)*resolution;

n_values = (max_edge - min_edge)/resolution + 1;

if n_values > MAX_N
    %TODO: Make the error more explicit
    error('Given specified resolution there are too many data bins')
end

edges = min_edge:resolution:max_edge;
bins  = edges(1:end-1) + resolution/2;

end