function initObject(obj,data)
%
%   seg_worm.stats.hist.initObject()
%
%   Old Name:
%   - seg_worm.util.histogram
%   - seg_worm.util.nanHistogram
%
%   This is essentially the constructor code. I moved it in here to avoid
%   the indenting.


%Remove NaN to speed up stats later ...

data(isinf(data) | isnan(data)) = [];

n_samples = length(data);

%NOTE: For empty sets, we use the default values ...
if n_samples == 0
    return
end

obj.samples = n_samples;

[obj.bins,edges] = h_computeBinInfo(data,obj.resolution);

% Compute the histogram counts for all the data.
%-------------------------------------------------------------
counts = histc(data, edges);
counts(end) = []; %Remove the extra bin at the end (for overflow)

obj.counts = counts;
obj.pdf    = counts./sum(counts);
obj.mean   = mean(data);

if n_samples == 1
    obj.std = 0;
else
    %I couldn't resist optimizing this ...
    obj.std = sqrt(1/(n_samples-1)*sum((data - obj.mean).^2));
end

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

%If we have a singular value, then we will get a singular edge, which isn't
%good for binning. We always need to make sure that our data is bounded by
%a high and low end. Given how hist works (it is inclusive on the low end,
%when we only have one edge we add a second edge by increasing the high
%end, NOT by decreasing the low end.
%
%i.e. in Matlab you can't bound 3 by having edges at 2 & 3, the edges would
%need to be at 3 & 4

if min_edge == max_edge
   max_edge = min_edge + resolution; 
end

n_values = (max_edge - min_edge)/resolution + 1;

if n_values > MAX_N
    %TODO: Make the error more explicit
    error('Given specified resolution there are too many data bins')
end

edges = min_edge:resolution:max_edge;
bins  = edges(1:end-1) + resolution/2;

end