classdef hist < handle
    %
    %   Class:
    %   seg_worm.stats.hist
    %
    %   Old Names: (NOTE: The last part of the names are the SegWorm names)
    %   - seg_worm.w.stats.worm2histogram 
    %   - seg_worm.util.histogram
    %   - seg_worm.util.nanHistogram
    %   - seg_worm.w.stats.addWormHistograms
    %
    %   Current Status:
    %   -----------------------------------------------------------------
    %   Working but could be optimized a bit and there are some missing
    %   features. At this point however it will work for computing the
    %   statistics objects. Documentation could also be improved.
    %
    %
    %   Questions:
    %   -----------------------------------------------------------------
    %
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
    %   File Dependencies
    %   ------------------------------------------------------------------
    %
    %
    %
    
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
        units               %units associated with the bin values
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
        
        resolution      %Bin resolution
        is_zero_bin     %This might be useless ...
        is_signed       %
        
        %--------------------------------------------------------------
        
        pdf     = NaN %[1 x n_bins] %Probability density value for each bin
        bins    = NaN %[1 x n_bins] %The center of each bin
        
        
        counts    = 0 %[n_videos x bins] %The # of values in each bin
        n_samples = 0 %[n_videos x 1] # of samples for each video, TODO: This
        %needs to be clarified, I also don't like the name ...
        
        
        mean_per_video = NaN  %[n_videos x 1]
        std_per_video  = NaN  %[n_videos x 1]
        
        
        %{
        Not included yet ...
        p_normal, only for n_valid_measurements >= 3
        [~,cur_s.p_normal]  = seg_worm.fex.swtest(cur_h_e.mean, 0.05, 0);
        
        q_normal - 
        %}
        
    end
    
    properties (Dependent)
        mean 
        std  %Standard deviation of means
        n_valid_measurements
    end
    
    properties (Hidden, Dependent)
        n_videos   %# of videos that the object contains ...
        first_bin
        last_bin
        n_bins
    end
    
    methods
        function value = get.mean(obj)
           value = nanmean(obj.mean_per_video); 
        end
        function value = get.std(obj)
           value = nanstd(obj.mean_per_video); 
        end
        function value = get.n_valid_measurements(obj)
           value = sum(~isnan(obj.mean)); 
        end
        function value = get.n_videos(obj)
            value = length(obj.mean);
        end
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
        function objs = hist(feature_file_paths)
            %
            %
            %   obj = seg_worm.stats.hist;
            %
            %   TODO: Allow loading from saved file as well.
            
            %Added this to allow deep copy code to work
            if nargin == 0
                return
            end
            
            objs = seg_worm.stats.hist.initObjects(feature_file_paths);
            
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
        all_hists = initObjects(feature_file_path)
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
            %This simplifies the merging process a bit. This is accomplished by
            %always setting bin edges at multiples of the resolution. This was
            %not done previously ...
            
            all_objs = [hist_cell_array{:}];
            
            n_videos_per_object = [all_objs(1,:).n_videos];
            
            if any(n_videos_per_object ~= 1)
                error('Multiple videos per object not yet implemented')
            end
            
            n_videos   = size(all_objs,2);
            n_features = size(all_objs,1);
            
            temp_results = cell(1,n_features);
            
            for iFeature = 1:n_features
                
                cur_feature_array = all_objs(iFeature,:);
                
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
                final_obj.bins      = new_bins;
                final_obj.counts    = new_counts;
                final_obj.n_samples       = cat(1,cur_feature_array.n_samples);
                final_obj.mean_per_video  = cat(1,cur_feature_array.std_per_video);
                final_obj.std_per_video   = cat(1,cur_feature_array.std_per_video);
                final_obj.pdf       = sum(final_obj.counts,2)./sum(final_obj.n_samples);
                
                %Hold onto final object for output
                temp_results{iFeature} = final_obj;
                
            end
            
            objs = [temp_results{:}];
            
        end
    end
    
end



%Old plotting code ...
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
