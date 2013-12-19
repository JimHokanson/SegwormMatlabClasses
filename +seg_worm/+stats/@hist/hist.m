classdef hist < handle
    %
    %   Class:
    %   seg_worm.stats.hist
    
    %Output normally contains three things
    %
    %Which of these are important?????
    %
    %
    %seg_worm.util.histogram
    
    properties
        %Identification
        %----------------------------------
        name %a.b.c.d
        type
        sign %-1,0,1,2 (2 all)
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
        
        pdf
        bins

        counts
        samples
        
        
        mean_all
        std_all
        
        %TODO: These should probably be for signed data only ...
        mean_abs
        std_abs
        
        mean_pos
        std_pos
        
        mean_neg
        std_neg
    
    %For events:
%        data: 0.1467
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
        function obj = hist(data,resolution,is_zero_bin,is_signed)
            %
            %
            %   obj = seg_worm.stats.hist(data,resolution,is_zero_bin,is_signed);
            
            initObject(obj,data,resolution,is_zero_bin,is_signed)
            
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
