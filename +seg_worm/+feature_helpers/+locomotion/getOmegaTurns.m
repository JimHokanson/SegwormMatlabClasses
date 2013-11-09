function getOmegaTurns()
%
%   JAH: I'm still working on this function ...
%

%Apparently I need to go even further back ...
%frontBlock ????

lastBlock = midBlock;
clear('midBlock');
secondLastBlock = frontBlock;
clear('frontBlock');
lenMid = length(secondLastBlock{1});
lenEnd = length(lastBlock{1});

% Define the part of mainBlock to be saved in the result file
dataSize = [lenMid-featureWindow+1, featureWindow+lenEnd];
% Get frame class
frameClass = [secondLastBlock{1}(lenMid-featureWindow+1:end), lastBlock{1}];
stageFlag = frameClass == 'm';

% Find omega upsilon bends
angleArray         = [secondLastBlock{5}(:,lenMid-featureWindow+1:end),lastBlock{5}];





[omegaFramesBlock, upsilonFramesBlock] = omegaUpsilonDetectCurvature(angleArray, stageFlag);




omegaFramesBlock   = omegaFramesBlock(dataSize(1):dataSize(2))';
upsilonFramesBlock = upsilonFramesBlock(dataSize(1):dataSize(2))';

% save data
featureData.omegaFrames   = [featureData.omegaFrames, omegaFramesBlock];
featureData.upsilonFrames = [featureData.upsilonFrames, upsilonFramesBlock]; 

% Ommega turns based on direction
%--------------------------------------------------------------------------
% First we will clean up the angle array
% We need to normalize angles accounting for dropped frames 
maxFrames  = round(fps/2);
lastAngle  = tailToHeadDirection(1);
gapCounter = 0;
wormDataFixed = nan(size(tailToHeadDirection));
for i = 2:length(tailToHeadDirection)
    if ~isnan(tailToHeadDirection(i))
        wormDataFixed(i-1) = tailToHeadDirection(i)-lastAngle;
        gapCounter = 0;
        
        lastAngle = tailToHeadDirection(i);
    else
        gapCounter = gapCounter + 1;
    end
    
    if gapCounter > maxFrames
        lastAngle = NaN;
    end
end

positiveJumps = find(wormDataFixed > 180) + 1;
negativeJumps = find(wormDataFixed < -180) + 1;

% subtract 2pi from remainging data after positive jumps
for j = 1:length(positiveJumps)
    tailToHeadDirection(positiveJumps(j):end) = tailToHeadDirection(positiveJumps(j):end) - 2*180;
end

% add 2pi to remaining data after negative jumps
for j = 1:length(negativeJumps)
    tailToHeadDirection(negativeJumps(j):end) = tailToHeadDirection(negativeJumps(j):end) + 2*180;
end

% get long NaN stretches
n = isnan(tailToHeadDirection);
% save start and end indices for the stretches
[start1, end1] = regexp( char(n+'A'), 'B{120,}', 'start', 'end' );

% interpolate missing data
tailToHeadDirection(isnan(tailToHeadDirection)) = interp1(...
    find(~isnan(tailToHeadDirection)),...
    tailToHeadDirection(~isnan(tailToHeadDirection)),...
    find(isnan(tailToHeadDirection)),'linear');

% return long NaN stretches back to NaN
if ~isempty(start1) && ~isempty(end1)
    for kk=1:size(start1,2)
        tailToHeadDirection(start1(kk):end1(kk)) = NaN;
    end
end
% %debug
% plot(abs(wormData));
% ylim(gca,[0, 360]);
% xlim(gca, [1, length(wormData)]);
1;


% Computer angle difference
windowSize = round(fps/4);
newData = nan(1, length(tailToHeadDirection));
for j=windowSize+1:length(tailToHeadDirection)- windowSize;
    newData(j) = tailToHeadDirection(j+windowSize) - tailToHeadDirection(j-windowSize);
end
tailToHeadDirectionChange = abs(newData/(windowSize*2));

tailToHeadDirectionChangeOmega = tailToHeadDirectionChange>3;

% Interpolate body angle
% get long NaN stretches
bodyAngle = featureData.angleArray;

bodyAngle(isnan(bodyAngle)) = interp1(find(~isnan(bodyAngle)),...
    bodyAngle(~isnan(bodyAngle)), find(isnan(bodyAngle)),'linear', 'extrap');


% Filter out all omega bend flags that are shorter than fps/4 frames (1/4 of a second)
n = tailToHeadDirectionChangeOmega == 1; 
% save start and end indices for the stretches
[start1, end1] = regexp( char(n+'A'), strcat('B{',num2str(round(fps/4)),',}'), 'start', 'end' );
tailToHeadDirectionChangeOmega = zeros(size(tailToHeadDirectionChangeOmega));
% return long NaN stretches back to NaN
if ~isempty(start1) && ~isempty(end1)
    for kk=1:size(start1,2)
        if mean(bodyAngle(start1(kk):end1(kk))) > 0
            tailToHeadDirectionChangeOmega(start1(kk):end1(kk)) = 1;
        else
            tailToHeadDirectionChangeOmega(start1(kk):end1(kk)) = -1;
        end
    end
end

featureData.omegaFrames = featureData.omegaFrames | tailToHeadDirectionChangeOmega;

% Now restore the sign
n = featureData.omegaFrames == 1; 
% save start and end indices for the stretches
[start1, end1] = regexp( char(n+'A'), strcat('B{',num2str(round(fps/4)),',}'), 'start', 'end' );
featureData.omegaFrames = zeros(size(featureData.omegaFrames));
% return long NaN stretches back to NaN
if ~isempty(start1) && ~isempty(end1)
    for kk=1:size(start1,2)
        if mean(bodyAngle(start1(kk):end1(kk))) > 0
            featureData.omegaFrames(start1(kk):end1(kk)) = 1;
        else
            featureData.omegaFrames(start1(kk):end1(kk)) = -1;
        end
    end
end


omegaFrames   = featureData.omegaFrames;

% Compute the omega frames.
omegaFramesDorsal = findEvent(omegaFrames, 1, [], true);
omegaFramesVentral = findEvent(omegaFrames, [], -1, true);

% Unify the ventral and dorsal turns.
omegaFrames = cat(2, omegaFramesVentral, omegaFramesDorsal);
isOmegaVentral = [true(1, length(omegaFramesVentral)), ...
    false(1, length(omegaFramesDorsal))];
if ~isempty(omegaFramesVentral) && ~isempty(omegaFramesDorsal)
    [~, orderI]     = sort([omegaFrames.start]);
    omegaFrames     = omegaFrames(orderI);
    isOmegaVentral  = isOmegaVentral(orderI);
end

% Compute the omega statistics.
[omegaEventStats, omegaStats] = events2stats(omegaFrames, fps, ...
    distance, [], 'interDistance');

% Add the turns ventral/dorsal side.
omegaFrames = omegaEventStats;
if ~isempty(omegaFrames)
    omegaCells = squeeze(struct2cell(omegaFrames));
    omegaCells{end+1,1} = [];
    for i = 1:size(omegaCells, 2)
        omegaCells{end, i} = isOmegaVentral(i);
    end
    omegaFieldNames          = fieldnames(omegaFrames);
    omegaFieldNames{end + 1} = 'isVentral';
    omegaFrames              = cell2struct(omegaCells, omegaFieldNames, 1);
end

% Reorganize everything for the feature file.
omegaFrequency = [];
omegaTimeRatio = [];
if ~isempty(omegaStats)
    omegaFrequency = omegaStats.frequency;
    omegaTimeRatio = omegaStats.ratio.time;
end
omegas = struct( ...
    'frames', omegaFrames, ...
    'frequency', omegaFrequency, ...
    'timeRatio', omegaTimeRatio);