function initObject(obj,exp_hist,ctl_hist)
%worm2StatsInfo  Compute worm statistics information and save it to a file.
%
%   seg_worm.stats.initObject
%
%   See Also:
%   seg_worm.stats.helpers.swtest
%

ALPHA = 0.05; %Not really used since we don't examine H, just p
TAIL  = 0;

obj.field              = exp_hist.field;
obj.name               = exp_hist.name; 
obj.short_name         = exp_hist.short_name;
obj.units              = exp_hist.units;
obj.feature_category   = exp_hist.feature_category;
obj.hist_type          = exp_hist.hist_type;
obj.motion_type        = exp_hist.motion_type;
obj.data_type          = exp_hist.data_type;

%zscore
if isnan(exp_hist.mean)
    obj.z_score_experiment = -Inf;
elseif isnan(ctl_hist.mean)
    obj.z_score_experiment = Inf;
else
    %This might need to be means_per_video, not the mean ...
    obj.z_score_experiment = (exp_hist.mean - ctl_hist.mean)/ctl_hist.std;
end

%TODO: Move this to the histogram, not here ...
%------------------------------------------------------------------------
p_fields  = {'p_normal_experiment' 'p_normal_control'};
hist_objs = {exp_hist ctl_hist};

for iObj = 1:2
    cur_field    = p_fields{iObj};
    cur_hist_obj = hist_objs{iObj};
    if cur_hist_obj.n_valid_measurements < 3
        obj.(cur_field) = NaN;
    else
        obj.(cur_field) = seg_worm.stats.helpers.swtest(cur_hist_obj.mean_per_video,ALPHA,TAIL);
    end
end

if isnan(exp_hist.mean) || isnan(ctl_hist.mean)
    %This is a literal trasnlation of the code (I think)
    %I'm a bit confused by it ...
    n_expt  = exp_hist.n_videos;
    n_total = n_expt + ctl_hist.n_videos;
    obj.p_w = seg_worm.stats.helpers.fexact(n_expt,n_total,n_expt,n_expt);
    obj.p_t = obj.p_w;
else
    obj.p_w = ranksum(exp_hist.valid_means,ctl_hist.valid_means);
end


%pWValues - these seem to be the real statistics used ...
%- exclusive - fexact  seg_worm.stats.helpers.fexact
%- ranksum






%Not reimplemented ...
%{
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
%}



%s - set zscore to 0
%in control but not experiment : -Inf
%in experiment but not control : Inf
%???? - What about if in neither?????


%TODO: Just pass specs object


%{

%This could be incorrect if we have control data as well ...
%------------------------------------------------------------
p_values_all = [stats_objs.p_normal];
valid_p_mask = ~isnan(p_values_all);

valid_p_values = p_values_all(valid_p_mask);

%:/ This is a call to the bioinformatics toolbox ...
[~,q_normal] = mafdr(valid_p_values);

q_normal_cell = num2cell(q_normal);

[stats_objs(valid_p_mask).q_normal] = deal(q_normal_cell{:});
%------------------------------------------------------------


%Things that are missing:
%- normalization if controls are present
%- q_values, for when controls are present and when they are not (above code)
%- significance code ...
%
%   What does seg_worm.w.stats.wormStats2Matrix do?
%       - puts everything in a matrix
%

%JAH: At this point ...

keyboard






zScoreMode = 'os';

% Are we permuting the data?
permutations = [];

% Are we displaying the progress?
isVerbose = false;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------






worm_field    = 'histogram';
control_field = 'histogram';



% Initialize the feature information.
dataInfo = seg_worm.w.stats.wormStatsInfo();

% Organize the worm features.
wormData    = h__worm2data(dataInfo, useWorm,    wormSortI,    wormData,    wormField,    isVerbose);
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







%% Normalize the data.
function wormData = h__normalizeData(wormData, controlData, isExclusive, zScoreMode)
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




%}


