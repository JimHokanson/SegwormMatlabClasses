function initObject(obj,data,resolution,is_zero_bin,is_signed)
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

if n_samples == 0
    return
end

[bins,edges] = h_computeBinInfo(data,resolution);

% Compute the histogram counts for all the data.
counts = histc(data, edges);
counts(end) = []; %Remove the extra bin at the end (for overflow)

pdfs = counts./sum(counts);

%JAH TODO: At this point ...


% Compute the means of the data sets.
means   = cellfun(@nanmean, data);
stdDevs = cellfun(@nanstd, data);
if isSigned
    
    % Compute the absolute value statisitics.
    absData    = cellfun(@(x) abs(x), data, 'UniformOutput', false);
    absMeans   = cellfun(@(x) nanmean(x), absData);
    absStdDevs = cellfun(@(x) nanstd(x), absData);
    
    % Compute the positive value statisitics.
    posData    = cellfun(@(x) x(x > 0), data, 'UniformOutput', false);
    posMeans   = cellfun(@(x) nanmean(x), posData);
    posStdDevs = cellfun(@(x) nanstd(x), posData);
    
    % Compute the negative value statisitics.
    negData    = cellfun(@(x) x(x < 0), data, 'UniformOutput', false);
    negMeans   = cellfun(@(x) nanmean(x), negData);
    negStdDevs = cellfun(@(x) nanstd(x), negData);
end




%Data - What we care about ...
%===========================================

% Organize the histogram data sets.
histData.data.counts          = counts;
histData.data.samples(:,1)    = n_samples; %Yikes, I think is is to force a colum vector ...
histData.data.mean.all(:,1)   = means;
histData.data.stdDev.all(:,1) = stdDevs;
if isSigned
    
    % Compute the absolute value statisitics.
    histData.data.mean.abs(:,1)   = absMeans;
    histData.data.stdDev.abs(:,1) = absStdDevs;
    
    % Compute the positive value statisitics.
    histData.data.mean.pos(:,1)   = posMeans;
    histData.data.stdDev.pos(:,1) = posStdDevs;
    
    % Compute the negative value statisitics.
    histData.data.mean.neg(:,1)   = negMeans;
    histData.data.stdDev.neg(:,1) = negStdDevs;
end



% Organize the histogram.
histData.PDF  = pdfs;
histData.bins = bins;
histData.resolution = resolution;
histData.isZeroBin  = isZeroBin;
histData.isSigned   = isSigned;
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