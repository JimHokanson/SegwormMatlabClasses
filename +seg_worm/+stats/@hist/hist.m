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
    
    %{
    
    Comparison to old code:
    ------------------------------------------------------------------
    This version of the histogram differs significantly from the old
    SegWorm histogram code. Notably:
    
    1) There is much less reliance on saving values to disk. There is no
    need to actually save files to disk.
    2) Two types of data are missing from the histogram.
        - stats computed on all video data merged together
        - stats computed on the stats, i.e. the mean of the means
    3) No special code is given to "experimental" vs "control" data.
    4) Signed histograms are separate objects instead of different
    properties in one object. This leads to more hist objects for a given
    set of features, but makes the code more straighforward. 
    
    Functionally nothing should be lost by doing this.
    
    %}
    
    
    properties
        %Identification
        %--------------------------------------------------------------
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
        %--------------------------------------------------------------
        
        pdf     = NaN %[1 x n_bins] %probability density value for each bin
        bins    = NaN %[1 x n_bins] %the center of each bin

        
        counts  = 0 %[n_videos x bins] %The # of values in each bin
        samples = 0  %[n_videos x 1] # of samples for each video, TODO: This needs to be
        %clarified, I also don't like the name ...
        
        
        mean = NaN  %[n_videos x 1]
        std  = NaN  %[n_videos x 1]
    end
    
    properties (Hidden, Dependent)
       first_bin
       last_bin
       n_bins
    end
    
    methods
        function value = get.first_bin(obj)
           value = obj.bins(1);
        end
        function value = get.last_bin(obj)
           value = obj.bins(end);
        end
        function value = get.n_bins(obj)
           value = length(obj.bins);
        end
    end

    methods
        function obj = hist(data,specs,hist_type,motion_type,data_type)
            %
            %
            %   obj = seg_worm.stats.hist(data,resolution,is_signed);
            %
            %   Called by:
            %   seg_worm.stats.hist.createHistograms
            
            %Added this to allow deep copy code to work
            if nargin == 0
                return
            end
            
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
        function obj_out = createCopy(obj_in)
            %
            %   This is used to create a merged object without affecting 
            %   the originals ...
            %
            %   obj_out = seg_worm.stats.hist.createCopy(obj_in)
            %
            %   NOTE: The stats are not currently copied as this is
            %   primarily for merging objects where the stats will be
            %   overwritten.
            
            
            obj_out = seg_worm.stats.hist;
           
            obj_out.feature_category = obj_in.feature_category;
            obj_out.resolution       = obj_in.resolution;
            obj_out.is_zero_bin      = obj_in.is_zero_bin;
            obj_out.is_signed        = obj_in.is_signed;
            obj_out.name             = obj_in.name;
            obj_out.short_name       = obj_in.short_name;
            obj_out.units            = obj_in.units;
            
            obj_out.hist_type   = obj_in.hist_type;
            obj_out.motion_type = obj_in.motion_type;
            obj_out.data_type   = obj_in.data_type;
            
           
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
        %
        %   Inputs
        %   ==============================================================
        %
        %   Outputs
        %   ==============================================================
        %
        %
        %
        %
        %   Currently each object should only contain a single set of data
        %   (i.e. single video) prior to merging. This could be changed.
        

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
        %
        %edges 1: 1,2,3,4,5
        %edges 2: 3.5,4.5,5.5,6.5
        %
        %Instead we might have:
        %edges 2: 3,4,5,6,7,8
        %
        %This simplifies the merging process a bit
        
        %TODO: Add check for multiple videos in any object, this is not yet
        %supported ...
        
        n_videos   = length(hist_cell_array);
        n_features = length(hist_cell_array{1});
        
        temp_results = cell(1,n_features);
        
        for iFeature = 1:n_features 
        
        %Here we are indexing into arrays of objects that are encased in a
        %cell array:
        temp = cellfun(@(x) x(iFeature),hist_cell_array,'un',0);
        cur_feature_array = [temp{:}]; %cellfun doesn't support arrays of objects
                
        %Create an output object with same meta properties
        final_obj   =  cur_feature_array(1).createCopy();

        
        %Align all bins
        %----------------------------------------------------------------
        n_bins     = [cur_feature_array.n_bins];
        start_bins = [cur_feature_array.first_bin];
        min_bin    = min(start_bins);
        max_bin    = max([cur_feature_array.last_bin]);
        
        cur_resolution = final_obj.resolution;
        new_bins       = min_bin:cur_resolution:max_bin;
        
        %Colon operator was giving warnings about non-integer indices :/
        start_indices = round((start_bins - min_bin)./cur_resolution + 1);
        end_indices   = start_indices + n_bins - 1;
        
        new_counts = zeros(length(new_bins),n_videos);
        
        for iVideo = 1:n_videos
           cur_start = start_indices(iVideo);
           if ~isnan(cur_start)
           cur_end   = end_indices(iVideo);
           new_counts(cur_start:cur_end,iVideo) = cur_feature_array(iVideo).counts;
           end
        end
        
        %Update final properties
        %------------------------------------------------------------------
        final_obj.bins    = new_bins;
        final_obj.counts  = new_counts;
        final_obj.samples = [cur_feature_array.samples]';
        final_obj.mean    = [cur_feature_array.mean]';
        final_obj.std     = [cur_feature_array.std]';
        final_obj.pdf     = sum(final_obj.counts,2)./sum(final_obj.samples);
        
        %Hold onto final object for output
        temp_results{iFeature} = final_obj;
        
        end
        
        objs = [temp_results{:}];
                    
        end
    end
    
end

% % %% Normalize binned data.
% % function [normData,normBins] = normBinData(data, bins, resolution)
% % 
% % % Initialize the normalized data and bins.
% % normData = [];
% % normBins = [];
% % if isempty(data)
% %     return;
% % end
% % 
% % % Compute the normalized bins.
% % minBin = min(cellfun(@(x) x(1), bins));
% % if isnan(minBin)
% %     return;
% % end
% % maxBin   = max(cellfun(@(x) x(end), bins));
% % numBins  = round((maxBin - minBin) / resolution) + 1;
% % normBins = linspace(minBin, maxBin, numBins);
% % 
% % % Normalize the data.
% % normData = zeros(length(data), numBins);
% % for i = 1:length(data)
% %     
% %     % No data.
% %     if isnan(bins{i})
% %         continue;
% %     end
% %     
% %     % Copy the data.
% %     startI = round((bins{i}(1) - minBin) / resolution) + 1;
% %     endI   = round((bins{i}(end) - minBin) / resolution) + 1;
% %     normData(i, startI:endI) = data{i};
% % end
% % end



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
