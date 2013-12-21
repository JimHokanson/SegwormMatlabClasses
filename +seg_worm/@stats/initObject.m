function initObject(hist_files, control_files, varargin)
%worm2StatsInfo  Compute worm statistics information and save it to a file.
%
%   This appears to be a top level function.
%
%   Old Code: seg_worm.w.stats.worm2StatsInfo
%
%   seg_worm.w.stats.worm2StatsInfo(FILENAME, 
%                  WORMFILES, 
%                  WORMINFOFILTER, 
%                  WORMFEATFILTER,
%                  CONTROLFILES, 
%                  CONTROLINFOFILTER, 
%                  CONTROLFEATFILTER,
%                  ZSCOREMODE, 
%                  PERMUTATIONS, 
%                  ISVERBOSE)
%
%   TOOLBOX FUNCTIONS
%   -----------------------------------------------------------------------
%   mafdr - bioinformatics toolbox 
%   
%   Implements: [2] Storey, J.D. (2002). A
%   direct approach to false discovery rates. Journal of the Royal
%   Statistical Society 64(3), 479–498.
%
%   
%
%
%   QUESTIONS
%   -----------------------------------------------------------------------
%   1) When only wormFiles are provided - no controls, what do the
%   statistical values represent?
%   2) When a file has controls, why are these ignored ????
%
%
%dataInfo = 
% 
% 726x1 struct array with fields:
% 
%     name
%     unit
%     title1
%     title2
%     title3
%     title1I
%     title2I
%     title3I
%     index
%     field
%     isMain
%     category
%     type
%     subType
%     sign
%
%wormData = 
% 
% 726x1 struct array with fields:
% 
%     zScore
%     dataMeans
%     dataStdDevs
%     dataSamples
%     mean
%     stdDev
%     samples
%     pNormal
%     qNormal
%
%
%   Outputs
%   =======================================================================
%       output_file_path - the file name for the worm statistics information;
%                  a file containing a structure with fields:
%
%                  wormInfo & controlInfo = the worm information
%
%                  *******************
%                  dataInfo: [726 x 1]
%
%                  name     = the feature's name
%                  units    = the feature's units
%                  title1   = the feature's 1st title
%                  title1I  = the feature's 1st title index
%                  title2   = the feature's 2nd title
%                  title2I  = the feature's 2nd title index
%                  title3   = the feature's 3rd title
%                  title3I   = the feature's 3rd title index
%                  field    = the feature's path; a struct where:
%
%                             histogram  = the histogram data path
%                             statistics = the statistics data path
%
%                  index    = the feature's field index
%                  isMain   = is this a main feature?
%                  category = the feature's category, where:
%
%                             m = morphology
%                             s = posture (shape)
%                             l = locomotion
%                             p = path
%
%                  type     = the feature's type, where:
%
%                             s = simple data
%                             m = motion data
%                             d = event summary data
%                             e = event data
%                             i = inter-event data
%
%                  subType  = the feature's sub-type, where:
%
%                             n = none
%                             f = forward motion data
%                             b = backward motion data
%                             p = paused data
%                             t = time data
%                             d = distance data
%                             h = frequency data (Hz)
%
%                  sign     = the feature's sign, where:
%
%                             s = signed data
%                             u = unsigned data
%                             a = the absolute value of the data
%                             p = the positive data
%                             n = the negative data
%
%                  *******************
%                  wormData & controlData:
%
%                  zScore      = the z-score per feature
%                                (normalized to the controls)
%                                Note 1: if no controls are present, the
%                                 zScore is left empty
%                                Note 2: if a feature is present in more
%                                 than one worm but absent in the controls,
%                                 the z-score is set to infinite;
%                                 conversely, if a feature is absent from
%                                 the worms but present in more than 1
%                                 control, the z-score is set to -infinite.
%                  dataMeans   = the feature data means
%                  dataStdDevs = the feature data standard deviations
%                  dataSamples = the feature data samples
%                  mean        = the mean per feature
%                  stdDev      = the standard deviation per feature
%                  samples     = the number of samples per feature
%                  pNormal     = the Shapiro-Wilk normality p-values
%                  qNormal     = the Shapiro-Wilk normality q-values
%
%
%                  significance.worm:
%
%                  Note: this field is only present with a control
%
%                  pValue              = the worm p-value
%                  qValue              = the worm q-value
%                  exclusiveFeaturesI  = the indices for exclusive features
%
%
%                  significance.features:
%
%                  Note: this field is only present with a control
%
%                  pTValue = the Student's t-test p-value(s), per feature
%                  qTValue = the Student's t-test q-value(s), per feature
%                  pWValue = the Wilcoxon rank-sum p-value(s), per feature
%                  qWValue = the Wilcoxon rank-sum q-value(s), per feature
%
%
%   Inputs
%   =======================================================================
%       wormFiles         - the worm histogram or statistics files
%       wormInfoFilter    - the worm information filtering criteria;
%                           a structure with any of the fields:
%
%              minFPS     = the minimum video frame rate (frames/seconds)
%              minTime    = the minimum video time (seconds)
%              maxTime    = the maximum video time (seconds)
%              minSegTime = the minimum time for segmented video (seconds)
%              minRatio   = the minimum ratio for segmented video frames
%              minDate    = the minimum date to use (DATENUM)
%              maxDate    = the maximum date to use (DATENUM)
%              years      = the years to use
%              months     = the months to use (1-12)
%              weeks      = the weeks to use (1-52)
%              days       = the days (of the week) to use (1-7)
%              hours      = the hours to use (1-24)
%              trackers   = the trackers to use (1-8)
%
%       wormFeatFilter    - the worm feature filtering criteria;
%                           a structure with the fields:
%
%               featuresI = the feature indices (see WORMDATAINFO)
%               minThr    = the minimum feature value (optional)
%               maxThr    = the maximum feature value (optional)
%               indices   = the sub indices for the features (optional)
%               subFields = the subFields for the features (optional)
%
%       controlFiles      - the control histogram or statistics files
%       controlInfoFilter - the control information filtering criteria
%       controlFeatFilter - the control feature filtering criteria
%       zScoreMode        - the z-score normalization mode; if one mode is
%                           provided, it's applied to both the worm and
%                           control; otherwise, the first mode applies to
%                           the worm and the second mode is applied to the
%                           control;
%                           the default is 'os'.
%
%                       o = normalize the worm to the Opposing group
%                       s = normalize the worm to itself (Same)
%                       p = normalize the worm to both groups (Population)
%
%       permutations     - the number of permutations to run for
%                          significance testing; if zero or empty, no
%                          permutations are run;
%                          the default is none
%       isVerbose        - verbose mode displays the progress;
%                          the default is yes (true)
%
%   See also:
%   seg_worm.w.stats.worm2histogram, 
%   seg_worm.w.stats.worm2stats, 
%   FILTERWORMINFO, 
%   seg_worm.w.stats.wormStatsInfo
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

zScoreMode = 'os';

% Are we permuting the data?
permutations = [];

% Are we displaying the progress?
isVerbose = false;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------



%JAH: Removed filters ...

worm_hist_data    = load();
control_hist_data = load();


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------






worm_field = 'histogram';
control_field = 'histogram';



% Initialize the feature information.
dataInfo = seg_worm.w.stats.wormStatsInfo();

% Organize the worm features.
wormData    = h__worm2data(dataInfo, useWorm, wormSortI, wormData, wormField, isVerbose);
controlData = h__worm2data(dataInfo, useControl, controlSortI, controlData, controlField, isVerbose);


% Normalize the worm features.
%--------------------------------------------------------------------------
if ~isempty(control_files)
    
    % Find exclusive features.
    isNaNWorm    = isnan([wormData.dataMeans]);
    isNaNControl = isnan([controlData.dataMeans]);
    isExclusive  = (all(isNaNWorm) & ~any(isNaNControl)) | (~any(isNaNWorm) & all(isNaNControl));
    exclusiveI   = find(isExclusive);
    
    % Normalize the worm features.
    wormData = h__normalizeData(wormData, controlData, isExclusive, zScoreMode(1));
    if length(zScoreMode) > 1
        controlData = h__normalizeData(controlData, wormData, isExclusive, zScoreMode(2));
    else
        controlData = h__normalizeData(controlData, wormData, isExclusive, zScoreMode(1));
    end
end

% Determine the significance of the worm features.
%--------------------------------------------------------------------------
significance.worm.pValue = [];
significance.worm.qValue = [];
significance.features = [];
if ~isempty(control_files)
    
    % Correct the Shapiro-Wilk for multiple testing.
    pNormal = [wormData.pNormal; controlData.pNormal];
    qNormal = nan(size(pNormal));
    nonNaNPNormal = pNormal(~isnan(pNormal));
    if ~isempty(nonNaNPNormal)
        qNormal(~isnan(pNormal)) = mafdr(nonNaNPNormal);
    end
    for i = 1:size(qNormal,2)
        wormData(i).qNormal    = qNormal(1,i);
        controlData(i).qNormal = qNormal(2,:);
    end
    
    % Are any of the features present in one strain but not the other?
    significance.worm.exclusiveFeaturesI = exclusiveI;
%     if ~isempty(significance.worm.exclusiveFeaturesI)
%         significance.worm.pValue = 0;
%         significance.worm.qValue = 0;
%     end
    
    % Compute the worm feature significance.
    pTValues = mattest([wormData.dataMeans]', [controlData.dataMeans]');
    pWValues = nan(length(dataInfo), 1);
    for i = 1:length(pWValues)
        
        % The features are exclusive, use Fisher's exact test.
        if isExclusive(i)
            numWorms = length(isNaNWorm(:,i));
            numTotal = numWorms + length(isNaNControl(:,i));
            p = fexact(numWorms, numTotal, numWorms, numWorms);
            pTValues(i) = p;
            pWValues(i) = p;
            
        % Use the Wilcoxon rank-sum test.
        elseif ~(all(isNaNWorm(:,i)) || all(isNaNControl(:,i)))
            wData = [wormData(i).dataMeans];
            cData = [controlData(i).dataMeans];
            pWValues(i) = ranksum(wData(~isNaNWorm(:,i)), ...
                cData(~isNaNControl(:,i)));
        end
    end
    
    % Correct for multiple testing.
    [~, qTValues] = mafdr(pTValues);
    [~, qWValues] = mafdr(pWValues);
    
    % Organize the worm feature significance.
    significance.features(length(dataInfo), 1).pTValue = [];
    significance.features(length(dataInfo), 1).pWValue = [];
    for i = 1:length(significance.features)
                    
        % The feature significance was measured.
        significance.features(i).pTValue = pTValues(i);
        significance.features(i).qTValue = qTValues(i);
        significance.features(i).pWValue = pWValues(i);
        significance.features(i).qWValue = qWValues(i);
    end
    
    % Compute the worm significance.
    significance.worm.pValue = min([significance.features.pWValue]);
    significance.worm.qValue = min([significance.features.qWValue]);
end

% Save the features.
if isempty(control_files)
    save(output_file_path, 'dataInfo', 'wormInfo', 'wormData', '-v7.3');
else
    save(output_file_path, 'dataInfo', 'wormInfo', 'wormData', 'controlInfo', 'controlData', 'significance', '-v7.3');
end
end



%% Organize the features.
function data = h__worm2data(info, useWorm, wormSortI, wormData, wormField)

%For the field of interest, grab:
%means
%std_devs
%samples
%
%Compute:
%mean
%stdev
%samples - this is the length of the means
%pNormal - seg_worm.fex.swtest - only do this when more than 3 samples ...
%qNormal - 



% Organize the worm features.
data(length(info),1).zScore      = [];
data(length(info),1).dataMeans   = [];
data(length(info),1).dataStdDevs = [];
data(length(info),1).dataSamples = [];
for i = 1:length(info)
    %each feature

    % The feature is not indexed.
    if isnan(info(i).index)
        for j = 1:length(wormData)
            
            % Get the feature data.
            %
            %   wtf = wormData{j}   
            %   
            %wtf.posture.wavelength.secondary.backward.statistics.data.mean.all
            %   -> [3x1 double]
            %
            %
            %   GETTING AN ERROR:
            %   Somewhere a NaN is not being inserted ...
            %
            dataMeans = seg_worm.util.getStructField(wormData{j}, info(i).field.(wormField{j}), true);
            
            % Organize the feature data.
            if isempty(dataMeans)
                data(i).dataMeans = cat(1, data(i).dataMeans, nan(sum(useWorm{j}), 1));
            else
                data(i).dataMeans = cat(1, data(i).dataMeans, dataMeans(useWorm{j}));
            end
            
            % Organize the feature data standard deviations and samples.
            if isempty(dataMeans) || info(i).type == 'd'
                data(i).dataStdDevs = cat(1, data(i).dataStdDevs, nan(sum(useWorm{j}), 1));
                data(i).dataSamples = cat(1, data(i).dataSamples, nan(sum(useWorm{j}), 1));
            else
                dataStdDevs = seg_worm.util.getStructField(wormData{j}, ...
                    h__field2StdDev(info(i).field.(wormField{j})), true);
                dataSamples = seg_worm.util.getStructField(wormData{j}, ...
                    h__field2Samples(info(i).field.(wormField{j})), true);
                data(i).dataStdDevs = ...
                    cat(1, data(i).dataStdDevs, dataStdDevs(useWorm{j}));
                data(i).dataSamples = ...
                    cat(1, data(i).dataSamples, dataSamples(useWorm{j}));
            end
        end
        
        % Sort the data.
        data(i).dataMeans   = data(i).dataMeans(wormSortI);
        data(i).dataStdDevs = data(i).dataStdDevs(wormSortI);
        data(i).dataSamples = data(i).dataSamples(wormSortI);
        
        % Compute the feature statistics.
        data(i).mean    = nanmean(data(i).dataMeans);
        data(i).stdDev  = nanstd(data(i).dataMeans);
        data(i).samples = sum(~isnan(data(i).dataMeans));
        if data(i).samples < 3
            data(i).pNormal = NaN;
        else
            [~, data(i).pNormal, ~] = seg_worm.fex.swtest(data(i).dataMeans, 0.05, 0);
        end
        
    % The feature is indexed.
    else
        for j = 1:length(wormData)
            
            
            %Basically all of this convoluted code obtains the field of
            %interest ...
            
            % Get the feature data.
            field = info(i).field.(wormField{j});
            indexI = strfind(field, ')');
            subField = [];
            eval(['subField = wormData{j}.' field(1:indexI) ';']);
            dataMeans = seg_worm.util.getStructField(subField, field((indexI + 2):end), true);
            
            % Organize the feature data.
            if isempty(dataMeans)
                data(i).dataMeans = ...
                    cat(1, data(i).dataMeans, nan(sum(useWorm{j}), 1));
            else
                data(i).dataMeans = ...
                    cat(1, data(i).dataMeans, dataMeans(useWorm{j}));
            end
            
            % Organize the feature data standard deviations and samples.
            if isempty(dataMeans) || info(i).type == 'd'
                data(i).dataStdDevs = cat(1, data(i).dataStdDevs, nan(sum(useWorm{j}), 1));
                data(i).dataSamples = cat(1, data(i).dataSamples, nan(sum(useWorm{j}), 1));
            else
                dataStdDevs = seg_worm.util.getStructField(subField, h__field2StdDev(field((indexI + 2):end)), true);
                dataSamples = seg_worm.util.getStructField(subField, h__field2Samples(field((indexI + 2):end)), true);
                data(i).dataStdDevs = cat(1, data(i).dataStdDevs, dataStdDevs(useWorm{j}));
                data(i).dataSamples = cat(1, data(i).dataSamples, dataSamples(useWorm{j}));
            end
        end
        
        % Sort the data.
        data(i).dataMeans   = data(i).dataMeans(wormSortI);
        data(i).dataStdDevs = data(i).dataStdDevs(wormSortI);
        data(i).dataSamples = data(i).dataSamples(wormSortI);

        % Compute the feature statistics.
        data(i).mean    = nanmean(data(i).dataMeans);
        data(i).stdDev  = nanstd(data(i).dataMeans);
        data(i).samples = sum(~isnan(data(i).dataMeans));
        if data(i).samples < 3
            data(i).pNormal = NaN;
        else
            [~, data(i).pNormal, ~] = seg_worm.fex.swtest(data(i).dataMeans, 0.05, 0);
        end
    end
end

% Correct for multiple testing.
pAllValues = [data.pNormal];
pValues = pAllValues(~isnan(pAllValues));
qNormal = [];
if ~isempty(pValues)
    %Bioinformatics toolbox :/
    %JAH TEMP FIX:
    %-------------------------------
    %???What are the defaults ????
    
    qNormal = pValues;
    %[~, qNormal] = mafdr(pValues);
end
qAllValues = nan(length(data), 1);
qAllValues(~isnan(pAllValues)) = qNormal;
for i = 1:length(data)
    data(i).qNormal = qAllValues(i);
end
end



%% Normalize the data.
function wormData = h__normalizeData(wormData, controlData, isExclusive, ...
    zScoreMode)
switch zScoreMode
    
    % Normalize the worm to the Opposing group.
    case 'o'
        for i = 1:length(wormData)
            
            % Does the worm have any exclusive features?
            wormData(i).zScore  = NaN;
            if isExclusive(i)
                
                % The worm has a feature not present in its control.
                if wormData(i).samples > controlData(i).samples
                    wormData(i).zScore = inf;
                
                % The worm lacks a feature present in its control.
                else
                    wormData(i).zScore = -inf;
                end
                
            % Normalize the worm to its control.
            elseif controlData(i).stdDev > 0
                wormData(i).zScore = (wormData(i).mean - controlData(i).mean) / controlData(i).stdDev;
            end
        end
        
    % Normalize the worm to itself (Same).
    case 's'
        for i = 1:length(wormData)
            wormData(i).zScore = 0;
        end
        
    % Normalize the worm to both groups (Population).
    case 'p'
        for i = 1:length(wormData)
            
            % Combine the data.
            allData = cat(1, wormData(i).dataMeans, controlData(i).dataMeans);
            
            % Normalize the worm.
            wormData(i).zScore  = NaN;
            zMean   = nanmean(allData);
            zStdDev = nanstd(allData);
            if zStdDev > 0
                wormData(i).zScore  = (wormData(i).mean - zMean) / zStdDev;
            end
        end
end
end



%% Convert a data mean field to its standard deviation equivalent.
function field = h__field2StdDev(field)

% Find the last "mean".
meanStr   = '.mean.';
stdDevStr = '.stdDev.';
i = strfind(field, meanStr);

% The "mean" cannot be found.
if isempty(i)
    field = [];
    
% Replace the mean with the standard deviation.
else
    field = [field(1:(i(end) - 1)) stdDevStr field((i(end) + length(meanStr)):end)];
end
end



%% Convert a data mean field to its samples equivalent.
function field = h__field2Samples(field)

% Find the last "mean".
meanStr    = '.mean.';
samplesStr = '.samples';
i = strfind(field, meanStr);

% The "mean" cannot be found.
if isempty(i)
    field = [];
    
% Replace the mean with the samples.
else
    field = [field(1:(i(end) - 1)) samplesStr];
end
end
