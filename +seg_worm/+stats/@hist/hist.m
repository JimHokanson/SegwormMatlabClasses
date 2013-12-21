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
    
    properties
        %Identification
        %----------------------------------
        name 
        short_name
        units
        feature_category %posture, movement, etc ???
        motion_type %'all' 'forward' 'paused' 'backward'
        data_type   %'all' 'absolute' 'positive' 'negative'
        resolution
        is_zero_bin
        is_signed
        
        %motion
        %--------------------------------
        %direction
        %all_directions
        %forward
        %backward
        %paused
        
        pdf      %[1 x n_bins] %probability density value for each bin
        bins     %[1 x n_bins] %the center of each bin

        counts   %[n_videos x bins] %The # of values in each bin
        samples = 0  %# of samples for each video, TODO: This needs to be
        %clarified, I also don't like the name ...
        
        
        mean %[n_videos x 1]
        std  %[n_videos x 1]
        
        
        
    
    %For events:
%        data: 0.1467  %[n_videos x 1]
%     samples: 1
%        mean: 0.1467
%      stdDev: 0
        
    end
    
    %Events
    %----------------------
    %{
    
    worm.locomotion.motion.backward
    
        frequency: [1x1 struct]
                       data: 0.1467
                    samples: 1
                       mean: 0.1467
                     stdDev: 0
            ratio: [1x1 struct]
                .time
                       data: 0.1467
                    samples: 1
                       mean: 0.1467
                     stdDev: 0
                .distance
                       data: 0.1467
                    samples: 1
                       mean: 0.1467
                     stdDev: 0
             time: [1x1 struct] - .histogram
         distance: [1x1 struct] - .histogram
        interTime: [1x1 struct] - .histogram
    interDistance: [1x1 struct] - .histogram
    
    worm.locomotion.turns.upsilons
            frequency: [1x1 struct]
            timeRatio: [1x1 struct]
                 time: [1x1 struct]
            interTime: [1x1 struct]
        interDistance: [1x1 struct]
    
    %}
    
    methods
        function obj = hist(data,specs,motion_type,data_type)
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
            
            obj.motion_type = motion_type;
            obj.data_type   = data_type;
            
            %TODO: Pass in spec, inherit specs from common spec
            
            initObject(obj,data)
            
        end
    end
    
    methods (Static)
        createHistograms(feature_file_path)
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
