classdef hist < handle
    %
    %   Class:
    %   seg_worm.stats.hist
    %
    %   See Also:
    %   seg_worm.util.histogram
    %   seg_worm.w.stats.addWormHistograms 
    %
    %
    %   Main Code:
    %   initObject
    %   
    %   Questions:
    %   -----------------------------------------------------------------
    %   When calculating the stats, how are counts and samples used????
    %
    %   TODO:
    %   --------------------
    %   cleanup property instantation and documentation ...
    %
    %   Missing Features:
    %   -----------------------------------------
    %   - saving to disk
    %   - version comparison
    %   - 
    %
    %
    
    properties
        %Identification
        %----------------------------------
        name 
        short_name
        units
        feature_category    %posture, locomotion, morphology, path
        
        hist_type  %'motion' 'simple' 'event'
        
        motion_type %'all' 'forward'    'paused'    'backward'
        %This is what the midbody of the worm is doing while these values
        %are obtained
        
        data_type   %'all' 'absolute'   'positive'  'negative'
        %This is an additional filter on the values of the data. Either all
        %values are included or:
        %
        %   - the absolute values are used
        %   - only positive values are used
        %   - only negative values are used
        
        resolution      %
        is_zero_bin     %This might be useless
        is_signed       %

        
        pdf      %[1 x n_bins] %probability density value for each bin
        bins     %[1 x n_bins] %the center of each bin

        
        counts   %[n_videos x bins] %The # of values in each bin
        samples = 0  %# of samples for each video, TODO: This needs to be
        %clarified, I also don't like the name ...
        
        
        mean = NaN  %[n_videos x 1]
        std  = NaN  %[n_videos x 1]
    end

    methods
        function obj = hist(data,specs,hist_type,motion_type,data_type)
            %
            %
            %   obj = seg_worm.stats.hist(data,resolution,is_signed);
            %
            %   Called by:
            %   seg_worm.stats.hist.createHistograms
            
            obj.feature_category = specs.feature_category;
            obj.resolution       = specs.resolution;
            obj.is_zero_bin      = specs.is_zero_bin;
            obj.is_signed        = specs.is_signed;
            obj.name             = specs.name;
            obj.short_name       = specs.short_name;
            obj.units            = specs.units;
            
            obj.hist_type   = hist_type;
            obj.motion_type = motion_type;
            obj.data_type   = data_type;
                        
            initObject(obj,data)
        end
%         function stats = getStats(obj)
%            %
%            %    Old Name:
%            %    worm2StatsInfo
%            %
%         end
    end
    
    methods (Static)
        all_hists = createHistograms(feature_file_path)
        function objs = mergeObjects(hist_cell_array)
        %
        %
        %   objs = seg_worm.stats.hist.mergeObjects(hist_cell_array)
        
        %??? - what gets merged, see addWormHistograms
        
        %Merging involves:
        %- adding samples
        %- checking hists for consistency of props
        %- concatenate props that are by n_videos
        %
        %??? what about pdf and bins
        
        %See end of class for how this works, needs to be rewritten ...
        
        %Eeek, this is for 708 features, not just 1!
        
        
        %The goal of this function is to go from n collections of 708
        %histogram summaries of features each, to one set of 708 histogram
        %summaries, that has n elements, one for each video
        %
        %i.e. from something like:
        %{a.b a.b a.b a.b} where .b may have size [1 x m]
        %
        %to:
        %
        %a.b, where .b is of size [n x m], in this case [4 x m]
        %
        %This requires merging histograms that are computed using different
        %bins. IMPORTANTLY, because of the way that we do the bin
        %definitions, bin edges will always match if they are close, i.e.
        %we'll NEVER have:
        %edges 1: 1,2,3,4,5
        %edges 2: 3.5,4.5,5.5,6.5
        %
        %Instead we might have:
        %edges 2: 3,4,5,6,7,8
        %
        %This simplifies the merging process a bit
        
        n_features = length(hist_cell_array{1});
        
        temp_results = cell(1,n_features);
        
        for iFeature = 1:n_features 
        
        temp = cellfun(@(x) x(iFeature),hist_cell_array,'un',0);
        cur_feature_array = [temp{:}]; %cellfun doesn't support arrays of objects
        
        keyboard
        
        %JAH: At this point ...
            
        all_bins   = cellfun(@(x) x.bins,hist_cell_array,'un',0);
        all_counts = cellfun(@(x) x.counts,hist_cell_array,'un',0);
        
        min_bin = min(cellfun(@(x) x(1),all_bins));
        max_bin = max(cellfun(@(x) x(end),all_bins));
        
        
        
        keyboard
        
        [counts,bins] = normBinData(counts, bins, resolution);

        % Combine the PDFs.
        %==============================
        % % % data.bins = bins;
        % % % data.PDF = zeros(1, length(bins));
        % % % numSets = 0;
        % % % for i = 1:length(data.data.samples)
        % % %     if data.data.samples(i) > 0
        % % %         data.PDF = data.PDF + data.data.counts(i,:) ./ data.data.samples(i);
        % % %         numSets = numSets + 1;
        % % %     end
        % % % end

        
        end
        

        
        keyboard
            
        end
    end
    
end

%% Normalize binned data.
function [normData,normBins] = normBinData(data, bins, resolution)

% Initialize the normalized data and bins.
normData = [];
normBins = [];
if isempty(data)
    return;
end

% Compute the normalized bins.
minBin = min(cellfun(@(x) x(1), bins));
if isnan(minBin)
    return;
end
maxBin   = max(cellfun(@(x) x(end), bins));
numBins  = round((maxBin - minBin) / resolution) + 1;
normBins = linspace(minBin, maxBin, numBins);

% Normalize the data.
normData = zeros(length(data), numBins);
for i = 1:length(data)
    
    % No data.
    if isnan(bins{i})
        continue;
    end
    
    % Copy the data.
    startI = round((bins{i}(1) - minBin) / resolution) + 1;
    endI   = round((bins{i}(end) - minBin) / resolution) + 1;
    normData(i, startI:endI) = data{i};
end
end



%{

 % Show the results.
            if verbose
                
                % Create a figure.
                h = figure;
                set(h, 'units', 'normalized', 'position', [0 0 1 1]);
                hold on;
                
                % Plot the histogram.
                if length(data) > 1
                    subplot(1, 2, 1);
                end
                histColor = [0 0 0];
                statColor = [1 0 0];
                plotHistogram(histData, 'Histogram', 'Value', 'Data', '', ...
                    histColor, statColor, 2);
                
                % Plot the contributing histograms.
                if length(data) > 1
                    
                    % Construct the contributing histograms.
                    histDatas(1:length(data)) = histData;
                    for i = 1:length(data)
                        
                        % Set the statistics.
                        histDatas(i).sets.samples    = 1;
                        histDatas(i).sets.mean       = histData.data.mean(i);
                        histDatas(i).sets.stdDev     = 0;
                        histDatas(i).data.counts     = histData.data.counts(i,:);
                        histDatas(i).data.samples    = histData.data.samples(i);
                        histDatas(i).data.mean       = histData.data.mean(i);
                        histDatas(i).data.stdDev     = histData.data.stdDev(i);
                        histDatas(i).allData.samples = histData.data.samples(i);
                        histDatas(i).allData.mean    = histData.data.mean(i);
                        histDatas(i).allData.stdDev  = histData.data.stdDev(i);
                        
                        % Set the PDF.
                        counts = histDatas(i).data.counts;
                        histDatas(i).PDF = counts ./ sum(counts);
                    end
                    
                    % Plot the contributing histograms.
                    subplot(1, 2, 2);
                    dataNames(1:length(data)) = {'Data'};
                    colors = lines(length(data));
                    seg_worm.util.plotHistogram(histDatas, 'Histogram', 'Value', dataNames, colors, colors);
                end
            endS

%}
