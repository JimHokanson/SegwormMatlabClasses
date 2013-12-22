% © Eviatar Yemini and Columbia University 2013
% You will not remove any copyright or other notices from the Software;
% you must reproduce all copyright notices and other proprietary notices
% on any copies of the Software.

function main()

%% WHAT SHOULD WE DO?

% Make a statistics spreadsheet (CSV format) for the experiments.
makeCSV = true;

% Make the TIF files to visualize the experiment results.
makeTIFs = true;

% Make a figure for each feature of every experiment.
makeFigures = true;


%% THE CODE BELOW IS FOR EXPERT USERS ONLY!

% Initialize the log file.
logFile = 'log.txt';

% Initialize the worm and feature filters.
% Note: these filters permit you to eliminate bad experiments with poor
% quality worm segmentations and/or improbable measurement results.
wormInfoFilter = [];
wormFeatureFilter = [];

% Initialize the quality of the TIF file.
% Note: higher quality trades off more time to completion and larger file
% sizes.
TIFquality = 2;

% Initialize the figures to make (by index). If empty, make all figures.
figureFeatureIndices = [];

% Initialize the control and experiment directory prefixes.
controlDirStr = 'Control-';
experimentDirStr = 'Experiment-';

% Initialize the directory names.
featuresDir = 'features'; % features directory
histogramsDir = 'histograms'; % histograms directory
statisticsDir = 'statistics'; % statistics directory
tifDir = 'TIFs'; % TIFs directory
figureDir = 'figures'; % figures directory

% Initialize the file names.
statisticsCSVFile = 'statistics.csv'; % statistics spreadsheet file
allStatisticsFile = 'all_statistics.mat'; % all statistics file
figureDataFile = 'figure_data.mat'; % figure data file

% Initialize the file suffixes.
featuresSuffix = '_features.mat'; % features suffix
histogramSuffix = '_histogram.mat'; % histogram suffix
histogramsSuffix = '_histograms.mat'; % all experiment histograms suffix
statisticsSuffix = '_statistics.mat'; % statistics suffix



%% PLEASE DO NOT EDIT BELOW THIS LINE!!!

%% SET THE WORKING DIRECTORY.
workingDir = uigetdir('~');
if ~workingDir
    return;
end
cd(workingDir);

%% LOG EVERYTHING.
if exist(['./' logFile], 'file')
    delete(logFile);
end
diary(logFile);

%% FIX THE RESOLUTION.
set(0, 'ScreenPixelsPerInch', 48);

%% IDENTIFY THE CONTROL AND EXPERIMENT(S)

% Identify the control strain and its directory.
dirs = dir([controlDirStr '*']);
if length(dirs) ~= 1
    error(['There must be one, and only one, directory of controls.' 10 ...
        'The format for the directory name is "Control-<STRAIN_NAME>"']);
end
controlDir = dirs.name;
controlName = strrep(controlDir, controlDirStr, '');
controlName = strtrim(controlName);
disp(['*** Control: ' controlName]);

% Identify the experiment strains and their directories.
dirs = dir([experimentDirStr '*']);
dirs = {dirs.name};
if length(dirs) < 1
    error(['There must be at least one directory of experiments.' 10 ...
        'The format for the directory names is "Experiment-<STRAIN_NAME>"']);
end
experimentDirs = cell(length(dirs),1);
experimentNames = cell(length(dirs),1);
for i = 1:length(dirs)
    experimentDirs{i} = dirs{i};
    experimentNames{i} = strrep(experimentDirs{i}, experimentDirStr, '');
    experimentNames{i} = strtrim(experimentNames{i});
    disp(['*** Experiment ' num2str(i) '/' num2str(length(dirs)) ' : ' ...
        experimentNames{i}]);
end



%% FIX THE STRAIN NAMES.

% Make the control feature files directory.
controlFeaturesDir = [controlDir '/' featuresDir];
if ~exist(['./' controlFeaturesDir], 'dir')
    mkdir(controlFeaturesDir);
end

% Move the control feature files to a subdirectory.
disp(['*** Moving the control feature files to "' controlFeaturesDir '".']);
featureFiles = dir([controlDir '/*' featuresSuffix]);
featureFiles = {featureFiles.name};
for i = 1:length(featureFiles)
    featureFile = featureFiles{i};
    srcFile = [controlDir '/' featureFile];
    dstFile = [controlFeaturesDir '/' featureFile];
    [success, msg, ~] = movefile(srcFile, dstFile);
    if ~success
        error(['Cannot move "' srcFile '" to "' dstFile '":' 10 msg]);
    end
end

% Fix the control's strain name in the features files.
disp(['*** Fixing the control strains to "' controlName '".']);
featureFiles = dir([controlFeaturesDir '/*' featuresSuffix]);
featureFiles = {featureFiles.name};
for i = 1:length(featureFiles)
    featureFile = [controlFeaturesDir '/' featureFiles{i}];
    load(featureFile, 'info');
    info.experiment.worm.strain = controlName;
    info.experiment.worm.genotype = controlName;
    save(featureFile, 'info', '-append');
end

% Move the experiment feature files to a subdirectory.
for i = 1:length(experimentDirs)
    
    % Make the experiment feature files directory.
    experimentDir = experimentDirs{i};
    experimentFeaturesDir = [experimentDir '/' featuresDir];
    if ~exist(['./' experimentFeaturesDir], 'dir')
        mkdir(experimentFeaturesDir);
    end
    
    % Move the experiment feature files to a subdirectory.
    disp(['*** Moving the "' experimentDir '" feature files to "' ...
        experimentFeaturesDir '".']);
    featureFiles = dir([experimentDir '/*' featuresSuffix]);
    featureFiles = {featureFiles.name};
    for j = 1:length(featureFiles)
        featureFile = featureFiles{j};
        srcFile = [experimentDir '/' featureFile];
        dstFile = [experimentFeaturesDir '/' featureFile];
        [success, msg, ~] = movefile(srcFile, dstFile);
        if ~success
            error(['Cannot move "' srcFile '" to "' dstFile '":' 10 msg]);
        end
    end

    % Fix the experiment's strain name in the features files.
    disp(['*** Fixing the "' experimentDir '" strains to "' ...
        experimentNames{i} '".']);
    featureFiles = dir([experimentFeaturesDir '/*' featuresSuffix]);
    featureFiles = {featureFiles.name};
    for j = 1:length(featureFiles)
        featureFile = [experimentFeaturesDir '/' featureFiles{j}];
        info = [];
        load(featureFile, 'info');
        info.experiment.worm.strain = experimentNames{i};
        info.experiment.worm.genotype = experimentNames{i};
        save(featureFile, 'info', '-append');
    end
end



%% CONVERT THE FEATURES TO HISTOGRAMS.

% Make the control histogram files directory.
controlHistogramsDir = [controlDir '/' histogramsDir];
if ~exist(['./' controlHistogramsDir], 'dir')
    mkdir(controlHistogramsDir);
end

% Convert the control feature files to histograms.
featureFiles = dir([controlFeaturesDir '/*' featuresSuffix]);
featureFiles = {featureFiles.name};
for i = 1: length(featureFiles)
    featureFile = [controlFeaturesDir '/' featureFiles{i}];
    histogramFile = [controlHistogramsDir '/' ...
        strrep(featureFiles{i}, featuresSuffix, histogramSuffix)];
    if ~exist(['./' histogramFile], 'file')
        disp(['*** Converting "' featureFile '" to "' histogramFile '".']);
        worm2histogram(histogramFile, featureFile, [], 1);
    end
end

% Aggregate the control histograms.
allControlHistogramsFile = [controlHistogramsDir '/' controlName ...
    histogramsSuffix];
if ~exist(['./' allControlHistogramsFile], 'file')
    disp(['*** Aggregating the "' controlName '" histograms into "' ...
        allControlHistogramsFile '".']);
    histogramFiles = dir([controlHistogramsDir '/*' histogramSuffix]);
    histogramFiles = {histogramFiles.name};
    for i = 1:length(histogramFiles)
        histogramFiles{i} = [controlHistogramsDir '/' histogramFiles{i}];
    end
    addWormHistograms(allControlHistogramsFile, histogramFiles, [], 0, 1);
end

% Convert the experiment feature files to histograms.
for i = 1:length(experimentDirs)
    
    % Make the experiment histogram files directory.
    experimentDir = experimentDirs{i};
    experimentHistogramsDir = [experimentDir '/' histogramsDir];
    if ~exist(['./' experimentHistogramsDir], 'dir')
        mkdir(experimentHistogramsDir);
    end
    
    % Convert the experiment feature files to histograms.
    experimentFeaturesDir = [experimentDir '/' featuresDir];
    featureFiles = dir([experimentFeaturesDir '/*' featuresSuffix]);
    featureFiles = {featureFiles.name};
    for j = 1: length(featureFiles)
        featureFile = [experimentFeaturesDir '/' featureFiles{j}];
        histogramFile = [experimentHistogramsDir '/' ...
            strrep(featureFiles{j}, featuresSuffix, histogramSuffix)];
        if ~exist(['./' histogramFile], 'file')
            disp(['*** Converting "' featureFile '" to "' ...
                histogramFile '".']);
            worm2histogram(histogramFile, featureFile, [], 1);
        end
    end
    
    % Aggregate the experiment histograms.
    allExperimentHistogramsFile = [experimentHistogramsDir '/' ...
        experimentNames{i} histogramsSuffix];
    if ~exist(['./' allExperimentHistogramsFile], 'file')
        disp(['*** Aggregating the "' experimentNames{i} ...
            '" histograms into "' allExperimentHistogramsFile '".']);
        histogramFiles = dir([experimentHistogramsDir '/*' ...
            histogramSuffix]);
        histogramFiles = {histogramFiles.name};
        for j = 1:length(histogramFiles)
            histogramFiles{j} = [experimentHistogramsDir '/' ...
                histogramFiles{j}];
        end
        addWormHistograms(allExperimentHistogramsFile, histogramFiles, ...
            [], 0, 1);
    end
end



%% CONVERT THE AGGREGATE HISTOGRAMS TO STATISTICS.

% Make the control statistics files directory.
controlStatisticsDir = [controlDir '/' statisticsDir];
if ~exist(['./' controlStatisticsDir], 'dir')
    mkdir(controlStatisticsDir);
end

% Convert the aggregate control histogram file to statistics.
controlStatisticsFile = [controlStatisticsDir '/' controlName ...
    statisticsSuffix];
if ~exist(['./' controlStatisticsFile], 'file')
    disp(['*** Converting "' allControlHistogramsFile '" to "' ...
        controlStatisticsFile '".']);
    worm2StatsInfo(controlStatisticsFile, allControlHistogramsFile, ...
        wormInfoFilter, wormFeatureFilter, [], [], [], [], [], 1);
end

% Convert the aggregate experiment histogram files to statistics.
for i = 1:length(experimentDirs)
    
    % Make the experiment statistics files directory.
    experimentDir = experimentDirs{i};
    experimentStatisticsDir = [experimentDir '/' statisticsDir];
    if ~exist(['./' experimentStatisticsDir], 'dir')
        mkdir(experimentStatisticsDir);
    end
    
    % Convert the aggregate experiment histogram files to statistics.
    experimentHistogramsDir = [experimentDir '/' histogramsDir];
    allExperimentHistogramsFile = [experimentHistogramsDir '/' ...
        experimentNames{i} histogramsSuffix];
    experimentStatisticsFile = [experimentStatisticsDir '/' ...
        experimentNames{i} statisticsSuffix];
    if ~exist(['./' experimentStatisticsFile], 'file')
        disp(['*** Converting "' allExperimentHistogramsFile '" to "' ...
            experimentStatisticsFile '".']);
        worm2StatsInfo(experimentStatisticsFile, ...
            allExperimentHistogramsFile, ...
            wormInfoFilter, wormFeatureFilter, ...
            allControlHistogramsFile, ...
            wormInfoFilter, wormFeatureFilter, [], [], 1);
    end
end

% Make the aggregate statistics files directory.
if ~exist(['./' statisticsDir], 'dir')
    mkdir(statisticsDir);
end

% Aggregate the statistics.
allStatisticsFile = [statisticsDir '/' allStatisticsFile];
if  ~exist(['./' allStatisticsFile], 'file')
    disp(['*** Aggregating all statistics into "' allStatisticsFile '".']);
    allStatisticsFiles = cell(length(experimentDirs) + 1, 1);
    allStatisticsFiles{1} = controlStatisticsFile;
    for i = 1:length(experimentDirs)
        allStatisticsFiles{i + 1} = [experimentDirs{i} '/' ...
            statisticsDir '/' experimentNames{i} statisticsSuffix];
    end
    wormStats2Matrix(allStatisticsFile, allStatisticsFiles, 1);
end



%% MAKE THE STATISTICS SPREADSHEET.

% Make the statistics spreadsheet.
if makeCSV
    
    % Make the statistics spreadsheet.
    statisticsCSVFile = [statisticsDir '/' statisticsCSVFile];
    if ~exist(['./' statisticsCSVFile], 'file')
        
        % Open the statistics spreadsheet.
        [fid, message] = fopen(statisticsCSVFile, 'w');
        if fid < 0
            error(['Cannot open "' statisticsCSVFile '" for writing: ' ...
                message]);
        end
        
        % Load the statistics.
        worm = [];
        load(allStatisticsFile, 'worm');
        
        % List all features.
        featureInfo = wormStatsInfo();
        
        % Initialize the printed data.
        % Print: feature number, feature name,
        % z-score, mean, standard deviation, SEM, sample size,
        % p(Shapiro Wilk), q(Shapiro Wilk)
        % p(T-test), q(T-test)
        % p(Wilcoxon Rank Sum), q(Wilcoxon Rank Sum)
        numFeaturePrintData = 2;
        extraFeatureColumns = 7;
        numStatisticsPrintData = 5 + 2 + 2 + 2;
        extraStatisticsColumns = 2;
        
        % Print space for the feature data.
        for i = 1:numFeaturePrintData
            fprintf(fid, ',');
        end
            
        % Print the extra space for the features.
        for i = 1:extraFeatureColumns
            fprintf(fid, ',');
        end
        
        % Print the strains header.
        strains = worm.info.strain;
        for i = 1:length(strains)
            
            % Print the strain.
            fprintf(fid, '"%s",', strains{i});
            
            % Print space for the statistics data.
            for j = 1:(numStatisticsPrintData - 1)
                fprintf(fid, ',');
            end
            
            % Print extra space for the statistics.
            for j = 1:extraStatisticsColumns
                fprintf(fid, ',');
            end
            
        end
        fprintf(fid, '\n');
        
        % Print the features header.
        fprintf(fid, '"Feature #",');
        fprintf(fid, '"Feature Data",');
        
        % Print the extra space for the features.
        for i = 1:extraFeatureColumns
            fprintf(fid, ',');
        end
        
        % Print the statistics header.
        for i = 1:length(strains)
            
            % Print the statistics information.
            fprintf(fid, '"Z-Score",');
            fprintf(fid, '"Mean",');
            fprintf(fid, '"S.D.",');
            fprintf(fid, '"SEM",');
            fprintf(fid, '"N",');
            fprintf(fid, '"p(Normal)",');
            fprintf(fid, '"q(Normal)",');
            fprintf(fid, '"p(T-Test)",');
            fprintf(fid, '"q(T-Test)",');
            fprintf(fid, '"p(Wilcoxon)",');
            fprintf(fid, '"q(Wilcoxon)",');

            % Print extra space for the statistics.
            for j = 1:extraStatisticsColumns
                fprintf(fid, ',');
            end
        end
        fprintf(fid, '\n');
        
        % Print the data.
        stats = worm.stats;
        sig = worm.sig;
        for i = 1:length(featureInfo)
            
            % Print the feature.
            feature = featureInfo(i);
            fprintf(fid, '%d,', i);
            fprintf(fid, '"%s",', [feature.name ' (' feature.unit ')']);
            
            % Print the extra space for the features.
            for j = 1:extraFeatureColumns
                fprintf(fid, ',');
            end
                
            % Print the strain data.
            for j = 1:length(strains)
                
                % Print the strain statistics.
                fprintf(fid, '%f,', stats.zScore(j,i));
                fprintf(fid, '%f,', stats.mean(j,i));
                fprintf(fid, '%f,', stats.stdDev(j,i));
                fprintf(fid, '%f,', ...
                    stats.stdDev(j,i) / sqrt(stats.samples(j,i)));
                fprintf(fid, '%d,', stats.samples(j,i));
                fprintf(fid, '%f,', stats.pNormal(j,i));
                fprintf(fid, '%f,', stats.qNormal.all(j,i));
                fprintf(fid, '%f,', sig.pTValue(j,i));
                fprintf(fid, '%f,', sig.qTValue.all(j,i));
                fprintf(fid, '%f,', sig.pWValue(j,i));
                fprintf(fid, '%f,', sig.qWValue.all(j,i));
                
                % Print extra space for the statistics.
                for k = 1:extraStatisticsColumns
                    fprintf(fid, ',');
                end
            end
            fprintf(fid, '\n');
        end
        
        % Close the file.
        fclose(fid);
    end
end



%% MAKE THE EXPERIMENT TIFS.

% Make the experiment tifs.
if makeTIFs
    
    % Make the TIF files directory.
    if ~exist(['./' tifDir], 'dir')
        mkdir(tifDir);
    end
    
    % Make the TIFs.
    for i = 1:length(experimentDirs)
        
        % Determine the TIF file name.
        experimentName = experimentNames{i};
        experimentDir = experimentDirs{i};
        TIFfile = [tifDir '/' experimentName '_vs_' controlName '.tif'];
        
        % Make the TIF file.
        if ~exist(['./' TIFfile], 'file')
            disp(['*** Creating "' TIFfile '".']);
            experimentHistogramsDir = [experimentDir '/' histogramsDir];
            allExperimentHistogramsFile = [experimentHistogramsDir '/' ...
                experimentName histogramsSuffix];
            experimentFeaturesDir = [experimentDir '/' featuresDir];
            worm2TIF(TIFfile, ...
                allExperimentHistogramsFile, experimentName, ...
                wormInfoFilter, wormFeatureFilter, ...
                experimentFeaturesDir, ...
                allControlHistogramsFile, controlName, ...
                wormInfoFilter, wormFeatureFilter, ...
                controlFeaturesDir, ...
                allStatisticsFile, false, TIFquality, true, true);
        end
    end
end



%% MAKE THE EXPERIMENT FIGURES.

% Make the experiment figures.
if makeFigures
    
    % Make the figure files directory.
    if ~exist(['./' figureDir], 'dir')
        mkdir(figureDir);
    end

    % Aggregate the figure-related data into a file.
    figureDataFile = [figureDir '/' figureDataFile];
    if ~exist(['./' figureDataFile], 'file')
        
        % Initialize the figure data.
        disp(['*** Aggregating figure data into "' figureDataFile '".']);
        figureData(length(experimentNames) + 1).name = [];
        figureData(length(experimentNames) + 1).data = [];
        figureData(length(experimentNames) + 1).sig = [];

        % Aggregate the control data.
        worm = [];
        load(allStatisticsFile, 'worm');
        wormData = [];
        load(controlStatisticsFile, 'wormData');
        figureData(1).name = controlName;
        figureData(1).data = wormData;
        figureData(1).sig.pN = worm.stats.pNormal(1,:);
        figureData(1).sig.qN = worm.stats.qNormal.all(1,:);
        figureData(1).sig.pT = [];
        figureData(1).sig.qT = [];
        figureData(1).sig.pW = [];
        figureData(1).sig.qW = [];
        
        % Aggregate the experiments data.
        for i = 1:length(experimentDirs)
            experimentName = experimentNames{i};
            experimentStatisticsDir = [experimentDirs{i} '/' statisticsDir];
            experimentStatisticsFile = [experimentStatisticsDir '/' ...
                experimentName statisticsSuffix];
            experimentI = find(cellfun(@(x) strcmp(x, experimentName), ...
                worm.info.strain), 1);
            wormData = [];
            load(experimentStatisticsFile, 'wormData');
            figureData(i + 1).name = experimentName;
            figureData(i + 1).data = wormData;
            figureData(i + 1).sig.pN = worm.stats.pNormal(experimentI,:);
            figureData(i + 1).sig.qN = worm.stats.qNormal.all(experimentI,:);
            figureData(i + 1).sig.pT = worm.sig.pTValue(experimentI,:);
            figureData(i + 1).sig.qT = worm.sig.qTValue.all(experimentI,:);
            figureData(i + 1).sig.pW = worm.sig.pWValue(experimentI,:);
            figureData(i + 1).sig.qW = worm.sig.qWValue.all(experimentI,:);
        end

        % Save the aggregate the figure-related data into a file.
        save(figureDataFile, 'figureData');
    end
    
    % Load the figure data.
    load(figureDataFile, 'figureData');
    
    % Initialize the figure labels.
    numLabels = length(experimentNames) + 1;
    tickSpace = 3;
    xTicks = 1:tickSpace:(numLabels * tickSpace);
    xTickLabels = cell(numLabels, 1);
    xTickLabels{1} = controlName;
    for i = 1:length(experimentNames)
        xTickLabels{i + 1} = experimentNames{i};
    end
    
    % Initialize the control's color scheme.
    controlStdDevColor = [1 .5 1] * .5;
    controlSEMColor = [.5 1 1] * .5;

    % List the features.
    featureInfo = wormStatsInfo();
    featureLabels = {featureInfo.name};
        
    % Create a list of all figures.
    if isempty(figureFeatureIndices)
        
        % List all features.
        figureFeatureIndices = 1:length(featureInfo);
        
        % Remove paused motion features indices.
        lowerLabels = lower(featureLabels);
        pausedI = cellfun(@(x) ~isempty(strfind(x, 'paused')), lowerLabels);
        crawlI = cellfun(@(x) ~isempty(strfind(x, 'crawl')), lowerLabels);
        figureFeatureIndices(pausedI & crawlI) = [];
    end

    % Make the figures.
    for i = 1:length(figureFeatureIndices)

        % Determine the figure file name.
        featureI = figureFeatureIndices(i);
        featureTitle = featureLabels{featureI};
        filename = featureTitle;
        filename = strrep(filename, '/', ' per ');
        filename = strrep(filename, '.', '');
        cutI = strfind(filename, '(');
        if ~isempty(cutI)
            filename = filename(1:(cutI - 2));
        end
        filename = [figureDir '/' num2str(featureI, '%03d') ' - ' ...
            strtrim(filename) '.eps'];

        % Make the figure.
        if ~exist(['./' filename], 'file')
            disp(['*** Creating figure ' num2str(i) '/' ...
                num2str(length(figureFeatureIndices)) ' "' ...
                featureTitle '".']);

            % Create the figure space.
            figureHandle = figure;
            set(figureHandle, 'units', 'centimeters', ...
                'PaperType', 'A4', 'position', [0 0 29.7 21], ...
                'Visible', 'off');
            set(figureHandle, 'units', 'normalized');
            hold on;
            
            % Plot the experiments.
            for j = 1:length(figureData)
                h = notBoxPlot(figureData(j).data(featureI).dataMeans, ...
                    xTicks(j));
                if j == 1
                    set([h.sdPtch], 'FaceColor', controlStdDevColor);
                    set([h.semPtch], 'FaceColor', controlSEMColor);
                end
            end
            
            % Add space for the statistics to the bottom of the figure.
            yLims = get(gca, 'YLim');
            yRange = diff(yLims);
            yMin = yLims(1) - yRange * 0.5;
            yMax = yLims(end);
            yStatistics = yLims(1) - yRange * 0.25;
            
            % Add the statistics to the bottom of the figure.
            for j = 1:length(figureData)
                data = figureData(j).data(featureI);
                statisticsLabel = cell(5,1);
                statisticsLabel{1} = ['N = ' num2str(data.samples)];
                statisticsLabel{2} = ['Mean = ' num2str(data.mean)];
                statisticsLabel{3} = ['S.D. = ' num2str(data.stdDev)];
                pN = figureData(j).sig.pN(featureI);
                qN = figureData(j).sig.qN(featureI);
                statisticsLabel{4} = ...
                    ['p(Normal) = ' num2str(pN) ' ' p2stars(pN)];
                statisticsLabel{5} = ...
                    ['q(Normal) = ' num2str(qN) ' ' p2stars(qN)];
                if j ~= 1
                    pT = figureData(j).sig.pT(featureI);
                    qT = figureData(j).sig.qT(featureI);
                    statisticsLabel{5} = ...
                        ['p(T-Test) = ' num2str(pT) ' ' p2stars(pT)];
                    statisticsLabel{6} = ...
                        ['q(T-Test) = ' num2str(qT) ' ' p2stars(qT)];
                    pW = figureData(j).sig.pW(featureI);
                    qW = figureData(j).sig.qW(featureI);
                    statisticsLabel{7} = ...
                        ['p(Wilcoxon) = ' num2str(pW) ' ' p2stars(pW)];
                    statisticsLabel{8} = ...
                        ['q(Wilcoxon) = ' num2str(qW) ' ' p2stars(qW)];
                end
                text(xTicks(j), yStatistics, statisticsLabel, ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle');
            end
            
            % Label the figure.
            title(upper(featureTitle));
            ylim([yMin yMax]);
            ylabel([featureTitle ' (' featureInfo(featureI).unit ')']);
            xlim([xTicks(1) - 1, xTicks(end) + 1]);
            set(gca, 'XTick', xTicks);
            set(gca, 'XTickLabel', xTickLabels);
            
            % Save the figure.
            print(figureHandle, '-cmyk', '-dpsc2', filename);
            close(figureHandle);
        end
    end
end



%% FINISH UP.
disp('*** All done!');
diary off;
end
