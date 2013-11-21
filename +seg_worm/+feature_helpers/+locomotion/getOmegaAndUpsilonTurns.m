function [omegas,upsilons] = getOmegaAndUpsilonTurns(nw,tail_to_head_direction)
%
%   JAH: I'm still working on this function ...
%
%   [omegas,upsilons] =
%   seg_worm.feature_helpers.locomotion.getOmegaAndUpsilonTurns(nw,tail_to_head_direction)


% %From supplemental of Nature Methods paper
% Turns. Omega and upsilon turn events are computed similarly to a previously 
% described method9 but using skeleton bends instead of a single head-midbody-tail 
% angle. Omega and upsilon turns are signed negatively whenever the worm?s 
% ventral side is sheltered within the concavity of its midbody bend. 
% 
% The worm bends (described in the section on ?Posture?) are used to find a 
% contiguous sequence of frames (interruptible by coiling and other segmentation 
% failures) wherein a large bend travels from the worm?s head, through its midbody, 
% to its tail. The worm?s body is separated into three equal parts from its head to its 
% tail. The mean supplementary angle is measured along each third. 
%
% For omega 
% turns, this angle must initially exceed 30° at the first but not the last third of the 
% body (the head but not the tail). 
%
% The middle third must then exceed 30°. 
%
% And finally, the last but not the first third of the body must exceed 30° (the tail but not 
% the head). 
%
% This sequence of a 30° mean supplementary angle, passing 
% continuously along the worm from head to tail, is labeled an omega turn event. 
%
% Upsilon turns are computed nearly identically but they capture all events that 
% escaped being labeled omega turns, wherein the mean supplementary angle 
% exceeded 15° on one side of the worm (the first or last third of the body) while not 
% exceeding 30° on the opposite end. 




FPS = 20; %TODO: Get this from a higher level ...

angle_array = nw.angles;
stage_flag  = nw.segmentation_status == 'm';

%???? How does this work ????
[omegaFrames, upsilonFrames] = seg_worm.feature_helpers.locomotion.omegaUpsilonDetectCurvature(angle_array, stage_flag);


tailToHeadDirectionChangeOmega = h__getDChangeOmega(FPS,tailToHeadDirection);



%I'm still working on rewriting/understanding this code ...


keyboard

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

%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================

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



%locomotion.turns.upsilons
%--------------------------------------------------------------------------
% Compute the upsilon frames.

%NOTE: This needs to be rewritten to use better event encapsulation ...

upsilonFramesDorsal  = findEvent(upsilonFrames, 1, [], true);
upsilonFramesVentral = findEvent(upsilonFrames, [], -1, true);

% Unify the ventral and dorsal turns.
%--------------------------------------------------------------------------
upsilonFrames = cat(2, upsilonFramesVentral, upsilonFramesDorsal);
isUpsilonVentral = [true(1, length(upsilonFramesVentral)), 
    false(1, length(upsilonFramesDorsal))];
if ~isempty(upsilonFramesVentral) && ~isempty(upsilonFramesDorsal)
    [~, orderI]      = sort([upsilonFrames.start]);
    upsilonFrames    = upsilonFrames(orderI);
    isUpsilonVentral = isUpsilonVentral(orderI);
end

% Compute the upsilon statistics.
%--------------------------------------------------------------------------
[upsilonEventStats, upsilonStats] = events2stats(upsilonFrames, fps, distance, [], 'interDistance');

% Add the turns ventral/dorsal side.
%--------------------------------------------------------------------------
upsilonFrames = upsilonEventStats;
if ~isempty(upsilonFrames)
    upsilonCells = squeeze(struct2cell(upsilonFrames));
    upsilonCells{end+1,1} = [];
    for i = 1:size(upsilonCells, 2)
        upsilonCells{end, i} = isUpsilonVentral(i);
    end
    upsilonFieldNames = fieldnames(upsilonFrames);
    upsilonFieldNames{end + 1} = 'isVentral';
    upsilonFrames = cell2struct(upsilonCells, upsilonFieldNames, 1);
end

upsilonFrequency = [];
upsilonTimeRatio = [];
if ~isempty(upsilonStats)
    upsilonFrequency = upsilonStats.frequency;
    upsilonTimeRatio = upsilonStats.ratio.time;
end
upsilons = struct( ...
    'frames', upsilonFrames, ...
    'frequency', upsilonFrequency, ...
    'timeRatio', upsilonTimeRatio);



end

function tailToHeadDirectionChangeOmega = h__getDChangeOmega(FPS,tailToHeadDirection)


%tail_to_head_direction

% Ommega turns based on direction
%--------------------------------------------------------------------------
% First we will clean up the angle array
% We need to normalize angles accounting for dropped frames 



%Some loop
%--------------------------------------------------------------------
max_frames = round(FPS/2);
is_good_th_direction_value = ~isnan(tailToHeadDirection);
lastAngle  = tailToHeadDirection(1);
gapCounter = 0;

wormDataFixed = nan(size(tailToHeadDirection));
for i = 2:length(tailToHeadDirection)    
    if is_good_th_direction_value(i)
        wormDataFixed(i-1) = tailToHeadDirection(i) - lastAngle;
        gapCounter = 0;
        lastAngle  = tailToHeadDirection(i);
    else
        gapCounter = gapCounter + 1;
    end
    if gapCounter > max_frames
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


% Computer angle difference
windowSize = round(fps/4);
newData = nan(1, length(tailToHeadDirection));
for j=windowSize+1:length(tailToHeadDirection)- windowSize;
    newData(j) = tailToHeadDirection(j+windowSize) - tailToHeadDirection(j-windowSize);
end

tailToHeadDirectionChange      = abs(newData/(windowSize*2));
tailToHeadDirectionChangeOmega = tailToHeadDirectionChange > 3;

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


end