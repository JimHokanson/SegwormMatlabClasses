%{

%????????????????????
%--------------------------------------------------------------------------
% Older version < 3 of segmentation had a bug in saving the failed frames.
% They have been indexed starting 0 not 1 (because of frame number being
% generated from time stamp rather than globalFrameCounter). To counter act
% it the indices for the frames that failed need to be added 1 to shift the
% failed frames by one and re-allign them. Here we will make a check for
% that and will raise a flag to add 1 in the upcoming loop
shiftFailedFrames = 0;

if ~isempty(failedFrames) && length(failedFrames(:,1)) > 2
    if sum(frameLabels(failedFrames(2:end,1))~='f') ~= 0
        shiftFailedFrames = 1;
    end
end
%--------------------------------------------------------------------------


%}