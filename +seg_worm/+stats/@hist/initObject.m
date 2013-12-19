function initObject(obj)
%
%   seg_worm.stats.hist.initObject()
%
%
%   This is essentially the constructor code.
%


%seg_worm.util.histogram
%seg_worm.util.nanHistogram

data(isinf(data)) = NaN;

n_samples = sum(~isnan(data));

if n_samples == 0
    return
end

%Compute the data range & padding
%------------------------------------------------
min_data = min(data);
max_data = max(data);

% Compute the padding.
if min_data < 0
    minPad = resolution - abs(rem(min_data, resolution));
else
    minPad = abs(rem(min_data, resolution));
end
if max_data < 0
    maxPad = abs(rem(max_data, resolution));
else
    maxPad = resolution - abs(rem(max_data, resolution));
end

% Translate the bins by half the resolution to create a zero bin.
% Note: we compute just enough padding to capture the data.
half_resolution = resolution / 2;
if isZeroBin
    if minPad > half_resolution
        minPad = minPad - half_resolution;
    else
        minPad = minPad + half_resolution;
    end
    if maxPad > half_resolution
        maxPad = maxPad - half_resolution;
    else
        maxPad = maxPad + half_resolution;
    end
end

% Compute the edge range.
minEdge = min_data - minPad;
maxEdge = max_data + maxPad;

% Compute the bins and their edges.
% Note: histc fills all bins with edges(k) <= data < edges(k + 1).
% The last bin is filled with data == edges(end).
% Therefore, we keep the last bin empty and throw it away to preserve
% equal bin sizes. For this reason the bin centers are spaced for
% their final composition (what they will look like after tossing away
% the empty last bin).
numBins = round((maxEdge - minEdge) / resolution);

%???? Why go in and not out????
%1 4
%3 - resolution 1
%
%4 - 1 => 3/1 => 3
%
%minEdge
bins    = linspace(minEdge + half_resolution, maxEdge - half_resolution, numBins);
edges   = bins - half_resolution;
edges(end + 1) = edges(end) + resolution;

% Fix the zero bin.
% Note: IEEE floating point issues may shift us just off zero.
if is_zero_bin
    [zeroBin, zeroI] = min(abs(bins));
    if zeroBin < half_resolution / 2
        bins(zeroI) = 0;
    end
end

% Compute the histogram counts for all the data.
counts = histc(data, edges);
if length(edges) > 1
    
    % Add the very last bin.
    if counts(1,end) > 0
        %???? When would this ever run????
        bins(end + 1) = bins(end) + resolution;
        
        % Strip off the empty last bin.
    else
        counts(end) = [];
    end
    
    % Strip off the empty first bin.
    if counts(1,1) == 0
        bins(1) = [];
        counts(1) = [];
    end
end

% Compute the counts for the data set.
allCounts(1,:) = counts;
if length(data) == 1
    pdfs(1,:) = counts ./ sum(counts);
    
    % Compute the normalized histogram for the data sets.
else
    counts = zeros(length(data), length(edges));
    pdfs   = zeros(length(data), length(edges));
    for i = 1:length(data)
        if ~empty(i)
            counts(i,:) = histc(data{i}, edges);
            pdfs(i,:)   = counts(i,:) ./ sum(counts(i,:));
        end
    end
    pdfs = mean(pdfs, 1);
    
    % Strip off the empty bin.
    if length(edges) > 1
        
        % Add a bin (this should never happen).
        if any(counts(:,end) > 0)
            bins(end + 1) = bins(end) + resolution;
            warning('histogram:LastBinNotEmpty', ...
                ['The last bin in the histogram is not empty like it ' ...
                'should be. Please contact the programmer)']);
            
            % Strip off the empty bin.
        else
            counts(:,end) = [];
            pdfs(:,end) = [];
        end
    end
end

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