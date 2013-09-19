function [frame_statuses,movesI,locations] = findStageMovement(obj)
%findStageMovement(infoFile, logFile, diffFile, verbose, guiHandle)
%FINDSTAGEMOVEMENT Find stage movements in a worm experiment.
%
%   [frames,movesI,locations] = findStageMovement(obj)
%
%   Algorithm:
%   See docs/Finding_Stage_Movement_Algorithm.m
%
%   [FRAMES INDICES LOCATIONS] = ...
%       FINDSTAGEMOVEMENT(INFOFILE, LOGFILE, DIFFFILE, VERBOSE, GUIHANDLE)
%
%   Input:
%       infoFile  - the XML file with the experiment information
%       logFile   - the CSV file with the stage locations
%       diffFile  - the MAT file with the video differentiation
%       verbose   - verbose mode 1 shows the results in a figure
%                   verbose mode 2 labels the stage movements in the figure
%       guiHandle - the GUI handle to use when showing the results;
%                   if empty, the results are shown in a new figure
%
%   Output:
%       frames    - a vector of frame status
%                   true  = the frame contains stage movement
%                   false = the frame does NOT contain stage movement
%                   NaN   = the original video frame was dropped
%                   Note: video frames are indexed from 0, Matlab indexes
%                   from 1, please adjust your calculations accordingly
%       movesI    - a 2-D matrix with, respectively, the start and end
%                   frame indices of stage movements
%       locations - the location of the stage after each stage movement
%
% See also:
%   VIDEO2DIFF
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

verbose = false;

%----------------------------------------------------------
% Is there a GUI handle?
if ~exist('guiHandle','var')
   guiHandle = []; 
end

fps        = obj.fps;
frameDiffs = obj.frame_diffs;

% Check the frame rate.
minFPS = .1;
maxFPS = 100;
if fps < minFPS || fps > maxFPS
    warning('video2Diff:WeirdFPS', [diffFile ' was recorded at ' ...
        num2str(fps) ' frames/second. An unusual frame rate']);
end

delayFrames = ceil(obj.info.delay_time * fps);
mediaTimes  = obj.media_times;
locations   = obj.locations;

% The media time must be initialized.
if isempty(mediaTimes) || mediaTimes(1) ~= 0
    error('findStageMovement:NoInitialMediaTime', ...
        'The first media time must be 0');
end

% If there's more than one initial media time, use the latest one.
%--------------------------------------------------------------------------
I = find(mediaTimes == 0,1,'last');
if I > 1
    % Save the spare 0 media time location in case the corresponding
    % stage-movement, frame difference occured after the video started.
    spareZeroTimeLocation = locations(I - 1,:);
    
    % Dump the extraneous 0 media times and locations.
    mediaTimes(1:(I - 1))  = [];
    locations(1:(I - 1),:) = [];
end

% No frame difference means the frame was dropped.
%???? - why, does this happen often?
frameDiffs(frameDiffs == 0) = NaN;

% Normalize the frame differences and shift them over one to align them
% with the video frames.
frameDiffs(2:(length(frameDiffs) + 1)) = frameDiffs / max(frameDiffs);
frameDiffs(1) = frameDiffs(2);

% Compute the global Otsu and small frame-difference thresholds.
% Note 1: we use the Otsu to locate frame-difference peaks corresponding to
% stage movement.
% Note 2: we use the small threshold to locate the boundaries between
% frame differences corresponding to stage movement and frame differences
% corresponding to a non-movement interval.
gOtsuThr    = graythresh(frameDiffs);
gSmallDiffs = frameDiffs(frameDiffs < gOtsuThr);
gSmallThr   = median(gSmallDiffs) + 3 * std(gSmallDiffs);


%--------------------------------------------------------------------------
% The log file doesn't contain any stage movements, quit early
if length(mediaTimes) < 2
    warning('findStageMovements:NoStageMovements', 'The stage never moves');
    
    % Are there any large frame-difference peaks?
    if gOtsuThr >= gSmallThr
        [~, indices] = seg_worm.util.maxPeaksDistHeight(frameDiffs, delayFrames, gOtsuThr);
        warning('findStageMovements:UnexpectedPeaks', ['There are ' ...
            num2str(length(indices)) ' large frame-difference peaks ' ...
            'even though the stage never moves']);
    end
    
    % Finish.
    frame_statuses = false(length(frameDiffs), 1);
    movesI = [0 0];
    return;
end




% Pre-allocate memory.
frame_statuses = false(length(frameDiffs), 1); % stage movement status for frames
movesI(1:length(mediaTimes), 1:2) = NaN; % stage movement indices
movesI(1,:) = 0;

%??????????????????
if verbose
    peaksI(1:length(mediaTimes)) = NaN; % stage movement frame peaks
    endPeaksI = []; % peaks after the last stage movement
    otsuThrs  = zeros(1, length(mediaTimes) + 1); % Otsu thresholds
    smallThrs = zeros(1, length(mediaTimes) - 1); % small thresholds
    timeOffs  = [mediaTimes, length(frameDiffs) / fps]; % offset media times
    
    % Have we matched all the media times to frame-difference stage movements?
    isMediaTimesMatched = true;
end

% Are there enough frames?
if sum(~isnan(frameDiffs)) < 2
    error('findStageMovement:IsufficientFrames', ...
        'The video must have at least 2, non-dropped frames');
end

% Compute the search boundary for the first frame-difference peak.
maxMoveFrames = delayFrames + 1;        % maximum frames a movement takes
maxMoveTime   = maxMoveFrames / fps;    % maximum time a movement takes
timeOff       = maxMoveTime;            % the current media time offset

peakI         = 1; % the current stage movement peak's index
prevPeakI     = 1; % the previous stage-movement peak's index
prevPeakEndI  = 1; % the previous stage-movement peak's end index

startI = 1; % the start index for our search
endI   = 2 * maxMoveFrames; % the end index for our search
if endI > length(frameDiffs)
    endI = length(frameDiffs);
end
searchDiffs = frameDiffs(startI:endI);

%--------------------------------------------------------------------------
% Is the Otsu threshold large enough?
otsuThr = graythresh(searchDiffs);
isOtsu = otsuThr > gOtsuThr; % false if no global Otsu
if ~isOtsu
    
    % Does the Otsu threshold separate the 99% of the small frame
    % differences from the large ones? And, if there is a global small
    % threshold, is the Otsu threshold larger?
    smallDiffs = searchDiffs(searchDiffs < otsuThr);
    isOtsu     = ~isempty(smallDiffs) && sum(~isnan(smallDiffs)) > 0 && ...
        (isnan(gSmallThr) || otsuThr > gSmallThr) && ...
        otsuThr >= median(smallDiffs) + 3 * std(smallDiffs);
    
    % Does the global Otsu threshold pull out any peaks?
    if ~isOtsu
        if ~isnan(gOtsuThr) && sum(searchDiffs > gOtsuThr) > 1
            otsuThr = gOtsuThr;
            isOtsu = true;
        end
    end
end

%--------------------------------------------------------------------------
% Are there any distinguishably large, frame-difference peaks?
if (verbose)
    peaksI(1) = peakI;
    otsuThrs(1) = gOtsuThr;
end

if isOtsu
    
    % Do the frame differences begin with a stage movement?
    indices = find(searchDiffs > otsuThr);
    firstPeakI = indices(1);
    if firstPeakI <= maxMoveFrames
        
        % Find the largest frame-difference peak.
        [~, peakI] = max(frameDiffs(1:maxMoveFrames));
        prevPeakI = peakI;
        
        % Compute the media time offset.
        timeOff = peakI / fps;
        if verbose
            peaksI(1)   = peakI;
            otsuThrs(1) = otsuThr;
            timeOffs(1) = timeOff;
        end
    end
    
    % Is there a still interval before the first stage movement?
    if peakI > 1
        i = peakI - 1;
        while i > 1
            if frameDiffs(i) < gSmallThr && frameDiffs(i - 1) < gSmallThr
                peakI = 1;
                break;
            end
            i = i - 1;
        end
    end
end


%Variables:
%---------------------------------------------
%peakI
%



%--------------------------------------------------------------------------
% We reached the end.
endI = peakI + maxMoveFrames;
if endI >= length(frameDiffs)
    prevPeakEndI = length(frameDiffs);
    
% Find a temporary front end for a potential initial stage movement.
else
    searchDiffs = frameDiffs(peakI:endI);
    
    % Does the search window contain multiple stage movements?
    if ~isnan(gOtsuThr) && ~isnan(gSmallThr)
        foundMove = false;
        for i = 1:length(searchDiffs)
            
            % We found a still interval.
            if ~foundMove && searchDiffs(i) < gSmallThr
                foundMove = true;
                
            % We found the end of the still interval, cut off the rest.
            elseif foundMove && searchDiffs(i) > gSmallThr
                searchDiffs = searchDiffs(1:(i - 1));
                break;
            end
        end
    end
    
    % Find a temporary front end for a potential initial stage movement.
    [minDiff,i]   = min(searchDiffs);
    peakFrontEndI = peakI + i - 1;
    
    % If the temporary front end's frame difference is small, try to push
    % the front end backwards (closer to the stage movement).
    if minDiff <= gSmallThr
        i = peakI;
        while i < peakFrontEndI
            if frameDiffs(i) <= gSmallThr
                peakFrontEndI = i;
                break;
            end
            i = i + 1;
        end
    
    % If the temporary front end's frame difference is large, try to
    % push the front end forwards (further from the stage movement).
    elseif minDiff >= gOtsuThr || (minDiff > gSmallThr && ...
            peakFrontEndI < endI && ...
            all(isnan(frameDiffs((peakFrontEndI + 1):endI))))
        peakFrontEndI = endI;
    end
    
    % Advance.
    prevPeakEndI = peakFrontEndI;
end


helper__superAwesome();




% Do the frame differences end with a stage movement?
if prevPeakEndI > length(frameDiffs)
    movesI(end,2) = length(frameDiffs);
    frame_statuses(movesI(end,1):end) = true;
    movesI(end + 1,:) = [length(frameDiffs) length(frameDiffs)] + 1;
    if verbose
        smallThrs(end + 1) = gSmallThr;
    end
    
% Find the front end for the last stage movement.
else
    
    % Is the Otsu threshold large enough?
    searchDiffs = frameDiffs(prevPeakEndI:end);
    otsuThr = graythresh(searchDiffs);
    isOtsu = otsuThr > gOtsuThr; % false if no global Otsu
    if ~isOtsu
        
        % Does the Otsu threshold separate the 99% of the small frame
        % differences from the large ones? And, if there is a global small
        % threshold, is the Otsu threshold larger?
        smallDiffs = searchDiffs(searchDiffs < otsuThr);
        isOtsu = ~isempty(smallDiffs) && sum(~isnan(smallDiffs)) > 0 && ...
            (isnan(gSmallThr) || otsuThr > gSmallThr) && ...
            otsuThr >= median(smallDiffs) + 3 * std(smallDiffs);
        
        % Does the global Otsu threshold pull out any peaks?
        if ~isOtsu
            if ~isnan(gOtsuThr) && sum(searchDiffs > gOtsuThr) > 1
                otsuThr = gOtsuThr;
                isOtsu = true;
            end
        end
    end
    
    % Are there any large frame difference past the last stage movement?
    isExtraPeaks = false;
    if ~isOtsu
        peakI = length(frameDiffs) + 1;
        peakBackEndI = length(frameDiffs);
        
    % There are too many large frame-difference peaks.
    else
        [~, indices] = maxPeaksDistHeight(searchDiffs, maxMoveFrames, otsuThr);
        isExtraPeaks = ~isempty(indices);
        if verbose
            endPeaksI = indices + prevPeakEndI - 1;
        end
        
        % Find the first large peak past the last stage movement.
        i = prevPeakEndI;
        while i < length(frameDiffs) && (isnan(frameDiffs(i)) || ...
                frameDiffs(i) < otsuThr)
            i = i + 1;
        end
        peakI = i;
        
        % Find a temporary back end for this large peak.
        % Note: this peak may serve as its own temporary back end.
        startI = max(peakI - maxMoveFrames, prevPeakEndI);
        [minDiff,i] = min(fliplr(frameDiffs(startI:peakI)));
        peakBackEndI = peakI - i + 1; % we flipped to choose the last min
        
        % If the temporary back end's frame difference is small, try to
        % push the back end forwards (closer to the stage movement).
        if minDiff <= prevSmallThr
            i = peakI - 1;
            while i > startI
                if frameDiffs(i) <= prevSmallThr
                    peakBackEndI = i;
                    break;
                end
                i = i - 1;
            end
            
        % If the temporary back end's frame difference is large, try to
        % push the back end backwards (further from the stage movement).
        elseif minDiff >= min(otsuThr, gOtsuThr) || ...
                (minDiff > gSmallThr && peakBackEndI > startI && ...
                all(isnan(frameDiffs(startI:(peakBackEndI - 1)))))
            peakBackEndI = startI;
        end
    end
            
    % Compute a threshold for stage movement.
    smallDiffs = frameDiffs(prevPeakEndI:peakBackEndI);
    smallThr = nanmean(smallDiffs) + 3 * nanstd(smallDiffs);
    if isnan(smallThr)
        smallThr = prevSmallThr;
    end
    if (verbose)
        smallThrs(end + 1) = smallThr;
    end
    
    % Find the front end for the last logged stage movement.
    i = prevPeakI;
    while i < peakI && ((isnan(frameDiffs(i)) || ...
            frameDiffs(i) > smallThr) && ...
            (isnan(frameDiffs(i + 1)) || ...
            frameDiffs(i + 1) > smallThr))
        i = i + 1;
    end
    movesI(end,2) = i - 1;
    prevPeakEndI = i - 1;
    
    % Mark the last logged stage movement.
    if size(movesI, 1) == 1
        frame_statuses(1:movesI(end,2)) = true;
    else
        frame_statuses(movesI(end,1):movesI(end,2)) = true;
    end
    
    % Are there any large frame-difference peaks after the last logged
    % stage movement?
    if isExtraPeaks
        warning('findStageMovement:TooManyPeaks', ...
            ['There are, approximately, ' num2str(length(indices)) ...
            ' large frame-difference peaks after the last stage' ...
            ' movement ends at ' num2str((movesI(end,2) - 1)/ fps, '%.3f') ...
            ' seconds (frame ' num2str(movesI(end,2) - 1) ')']);
    end
    
    % Find the back end for logged stage movements.
    i = peakI - 1;
    while i > prevPeakEndI && (isnan(frameDiffs(i)) || ...
            frameDiffs(i) > smallThr)
        i = i - 1;
    end
    movesI(end + 1,:) = [i + 1, length(frameDiffs) + 1];
    frame_statuses(movesI(end,1):end) = true;
end

% Are any of the stage movements considerably small or large?
if (~verbose || isMediaTimesMatched) && isExtraPeaks
    
    % Compute the stage movement sizes.
    movesI(i:end,:) = [];
    moveSizes = zeros(size(movesI, 1),1);
    for j = 2:(size(movesI, 1) - 1)
        
        moveDiffs = frameDiffs(movesI(j,1):movesI(j,2));
        moveSizes(j) = sum(moveDiffs(~isnan(moveDiffs)));
        
%         % Interpolate over NaNs.
%         moveDiffs = frameDiffs((movesI(j,1) - 1):(movesI(j,2) + 1));
%         moveDiffs(isnan(moveDiffs)) = ...
%             interp1(find(~isnan(moveDiffs)), ...
%             moveDiffs(~isnan(moveDiffs)), find(isnan(moveDiffs)), ...
%             'linear');
%         moveSizes(j) = sum(moveDiffs(~isnan(moveDiffs(2:end-1))));
    end
    
    % Compute the statistics for stage movement sizes.
    meanMoveSize = mean(moveSizes(2:end));
    stdMoveSize = std(moveSizes(2:end));
    smallMoveThr = meanMoveSize - 2.5 * stdMoveSize;
    largeMoveThr = meanMoveSize + 2.5 * stdMoveSize;
    
    % Are any of the stage movements considerably small or large?
    for i = 2:(size(movesI, 1) - 1)
        
        % Is the stage movement small?
        if moveSizes(i) < smallMoveThr
            
            % Report the warning.
            warning('findStageMovement:ShortMove', ...
                ['Stage movement ' num2str(i) ...
                ' at media time ' num2str(mediaTimes(i), '%.3f') ...
                ' seconds (frame ' ...
                num2str(round(mediaTimes(i) * fps)) ...
                '), spanning from ' ...
                num2str((movesI(i,1) - 1) / fps, '%.3f') ...
                ' seconds (frame ' num2str(movesI(i,1) - 1) ...
                ') to ' num2str((movesI(i,2) - 1) / fps, '%.3f') ...
                ' seconds (frame ' ...
                num2str(movesI(i,2) - 1) '), is considerably small']);
            
        % Is the stage movement large?
        elseif moveSizes(i) > largeMoveThr
            
            % Report the warning.
            warning('findStageMovement:LongMove', ...
                ['Stage movement ' num2str(i) ...
                ' at media time ' num2str(mediaTimes(i), '%.3f') ...
                ' seconds (frame ' ...
                num2str(round(mediaTimes(i) * fps)) ...
                '), spanning from ' ...
                num2str((movesI(i,1) - 1) / fps, '%.3f') ...
                ' seconds (frame ' num2str(movesI(i,1) - 1) ...
                ') to ' num2str((movesI(i,2) - 1) / fps, '%.3f') ...
                ' seconds (frame ' ...
                num2str(movesI(i,2) - 1) '), is considerably large']);
        end
    end
end

% Show the stage movements.

end

function helper__superAwesome()

keyboard

%--------------------------------------------------------------------------
% Match the media time-stage movements to the frame-difference peaks.
mediaTimeOff = 0; % the offset media time
prevOtsuThr  = gOtsuThr; % the previous small threshold
prevSmallThr = gSmallThr; % the previous small threshold
isShifted    = false; % have we shifted the data to try another alignment?
i = 1;
while i < length(mediaTimes)
    
    % Advance.
    i = i + 1;
    
    % Compute the offset media time.
    prevMediaTimeOff = mediaTimeOff;
    mediaTimeOff = mediaTimes(i) + timeOff;
    if verbose
        timeOffs(i) = mediaTimeOff;
    end
    
    % Compute the search boundary for matching frame-difference peaks.
    mediaTimeOffI = round(mediaTimeOff * fps);
    startI        = prevPeakEndI;
    endI          = max(startI + 2 * abs(mediaTimeOffI - startI), max(startI, mediaTimeOffI) + maxMoveFrames);
    if endI > length(frameDiffs)
        endI = length(frameDiffs);
    end
    searchDiffs = frameDiffs(startI:endI);
    
    % Is the Otsu threshold large enough?
    otsuThr = graythresh(searchDiffs);
    isOtsu  = otsuThr > prevSmallThr || otsuThr > gOtsuThr;
    if ~isOtsu
        
        % Does the Otsu threshold separate the 99% of the small frame
        % differences from the large ones?
        if isnan(prevSmallThr) || otsuThr > prevSmallThr || otsuThr > gSmallThr
            smallDiffs = searchDiffs(searchDiffs < otsuThr);
            isOtsu = ~isempty(smallDiffs) && ...
                sum(~isnan(smallDiffs)) > 0 && ...
                otsuThr >= median(smallDiffs) + 3 * std(smallDiffs);
        end
        
        % Try the global Otsu threshold or, if there is none, attempt to
        % use half the search window's maximum frame difference.
        if ~isOtsu
            
            % Try using half the search window's maximum frame difference.
            if isnan(gOtsuThr)
                otsuThr = max(searchDiffs) / 2;
                
                % Does the half-maximum threshold separate the 99% of the
                % small frame differences from the large ones?
                smallDiffs = searchDiffs(searchDiffs < otsuThr);
                isOtsu = ~isempty(smallDiffs) && ...
                    sum(~isnan(smallDiffs)) > 0 && ...
                    otsuThr >= median(smallDiffs) + 3 * std(smallDiffs);
                
            % Does the global Otsu threshold pull out any peaks?
            elseif sum(searchDiffs > gOtsuThr) > 0
                otsuThr = gOtsuThr;
                isOtsu = true;
                
            % Does the global Otsu threshold pull out any peaks?
            elseif sum(searchDiffs > prevOtsuThr) > 0
                otsuThr = prevOtsuThr;
                isOtsu = true;
            end
        end
    end
    if (verbose)
        otsuThrs(i) = otsuThr;
    end
    
    % If we're at the end, make sure we're using an appropriate threshold.
    if i == length(mediaTimes)
        
        % Does the threshold separate the 99% of the small frame
        % differences from the large ones?
        smallDiffs = searchDiffs(searchDiffs < otsuThr);
        isOtsu = ~isempty(smallDiffs) && sum(~isnan(smallDiffs)) > 0 && ...
            otsuThr >= median(smallDiffs) + 3 * std(smallDiffs);
    end
    
    % Match the media time stage movement to a peak.
    indices = [];
    if isOtsu
        
        % Compute and set the global thresholds.
        if isnan(gOtsuThr)
            
            % Use a small threshold at 99% of the small frame differences.
            smallDiffs = searchDiffs(searchDiffs < otsuThr);
            smallThr = median(smallDiffs) + 3 * std(smallDiffs);
            
            % Set the global thresholds.
            if otsuThr >= smallThr
                gOtsuThr = otsuThr;
                gSmallThr = smallThr;
                
                % Set the previous small threshold.
                if isnan(prevOtsuThr)
                    prevOtsuThr = otsuThr;
                    prevSmallThr = smallThr;
                end
                
            % Use the previous small threshold.
            elseif ~isnan(prevSmallThr)
                smallThr = prevSmallThr;
            end
            
        % Compute the local thresholds.
        else
            otsuThr = min(otsuThr, gOtsuThr);
            smallThr = max(prevSmallThr, gSmallThr);
            if smallThr > otsuThr
                smallThr = min(prevSmallThr, gSmallThr);
            end
        end
        
        % Does the search window contain multiple stage movements?
        foundMove = false;
        for j = 1:length(searchDiffs)
            
            % We found a stage movement.
            if ~foundMove && searchDiffs(j) > otsuThr
                foundMove = true;
                
            % We found the end of the stage movement, cut off the rest.
            elseif foundMove && searchDiffs(j) < smallThr
                searchDiffs = searchDiffs(1:(j - 1));
                break;
            end
        end
        
        % Find at least one distinguishably large peak.
        [~, indices] = ...
            maxPeaksDistHeight(searchDiffs, maxMoveFrames, otsuThr);
    end
    
    % We can't find any distinguishably large peaks.
    peakI = [];
    if isempty(indices)
        
        % Does the last stage movement occur after the video ends?
        if i == length(mediaTimes) && endI >= length(frameDiffs)
            
            % Does the last offset media time occur before the video ends?
            if mediaTimeOff < (length(frameDiffs) - 1) / fps
                warning('findStageMovement:LastPeak', ...
                    ['The search window for the last stage movement (' ...
                    num2str(i) ') at media time ' ...
                    num2str(mediaTimes(i), '%.3f') ...
                    ' seconds (frame ' ...
                    num2str(round(mediaTimes(i) * fps)) ...
                    ') offset to ' num2str(mediaTimeOff, '%.3f') ...
                    ' seconds (frame ' num2str(round(mediaTimeOff * fps)) ...
                    '), spanning from ' ...
                    num2str((startI - 1) / fps, '%.3f') ...
                    ' seconds (frame ' num2str(startI - 1) ...
                    ') to the last frame ' ...
                    num2str((endI - 1) / fps, '%.3f') ' seconds (frame ' ...
                    num2str(endI - 1) '), doesn''t have any'...
                    ' distinguishably large peaks. The peak probably' ...
                    ' occured after the video ended and, therefore,' ...
                    ' the last stage movement will be ignored.']);
            end
            
            % Ignore the last stage movement.
            mediaTimes(end)  = [];
            locations(end,:) = [];
            movesI(end,:)    = [];
            if verbose
                peaksI(end)    = [];
                otsuThrs(end)  = [];
                smallThrs(end) = [];
                timeOffs(end)  = [];
            end
            break;
        end
        
        % Report the warning.
        warning('findStageMovement:NoPeaks', ...
            ['The search window for stage movement ' num2str(i) ...
            ' at media time ' num2str(mediaTimes(i), '%.3f') ...
            ' seconds (frame ' num2str(round(mediaTimes(i) * fps)) ...
            ') offset to ' num2str(mediaTimeOff, '%.3f') ...
            ' seconds (frame ' num2str(round(mediaTimeOff * fps)) ...
            '), spanning from ' num2str((startI - 1) / fps, '%.3f') ...
            ' seconds (frame ' num2str(startI - 1) ') to ' ...
            num2str((endI - 1) / fps, '%.3f') ' seconds (frame ' ...
            num2str(endI - 1) '), doesn''t have any distinguishably' ...
            ' large peaks']);
        
    % Use the first peak.
    else
        peakI = indices(1) + startI - 1;
        if verbose
            peaksI(i) = peakI;
        end
        
        % Is the current offset media time further from the frame-
        % difference stage movement than the previous offset media time?
        peakTime = (peakI - 1) / fps;
        timeDiff = mediaTimeOff - peakTime;
        prevTimeDiff = prevMediaTimeOff - peakTime;
        if i > 2 && (abs(prevTimeDiff) > maxMoveTime || ...
                abs(timeDiff) > maxMoveTime) && ...
                mediaTimeOff > prevMediaTimeOff && ...
                abs(timeDiff / prevTimeDiff) > 2
            warning('findStageMovement:FarPeak', ...
                ['Stage movement ' num2str(i) ' (at media time ' ...
                num2str(mediaTimes(i), '%.3f') ' seconds) offset to ' ...
                num2str(mediaTimeOff, '%.3f') ...
                ' seconds, has its frame-difference peak at ' ...
                num2str(peakTime, '%.3f') ' seconds (frame ' ...
                num2str(peakI - 1) '), an error of ' ...
                num2str(timeDiff, '%.3f') ' seconds.' ...
                ' The previous media time, offset to ' ...
                num2str(prevMediaTimeOff, '%.3f') ' seconds, is closer' ...
                ' with an error of only ' num2str(prevTimeDiff, '%.3f') ...
                ' seconds (less than half the current media time''s' ...
                ' error). Therefore, we probably have either a false' ...
                ' peak, a shifted misalignment, or an abnormally long delay']);
            
            % Ignore this wrong peak.
            peakI = [];
        end
    end
    
    % Can we realign (shift) the stage movements and media times?
    if isempty(peakI)
        lastMoveTime = movesI(i - 1,1) / fps;
        isShiftable = true;
        if isShifted
            isShiftable = false;
            
        % Shift the media times forward.
        elseif i > 2 && abs(mediaTimes(i - 2) - lastMoveTime) < ...
                abs(mediaTimes(i) - lastMoveTime)
            
            % Would a time shift align the media times with the
            % frame-difference stage movements?
            for j = 2:(i - 2)
                
                % Compute the error from the predicted time.
                offset =  movesI(j,1) / fps - mediaTimes(j - 1);
                predictedTime = mediaTimes(j) + offset;
                moveTime =  movesI(j + 1,1) / fps;
                timeDiff = abs(predictedTime - moveTime);
                
                % Compute the interval between the media times.
                mediaDiff = mediaTimes(j) - mediaTimes(j - 1);
                
                % Is the error in the predicted time greater than
                % the interval between media times?
                if timeDiff > mediaDiff
                    isShiftable = false;
                    break;
                end
            end
            
            % Time cannot be shifted due to misalignment between the media
            % times and frame-difference stage movements.
            if ~isShiftable
                warning('findStageMovement:TimeShiftAlignment', ...
                    ['Time cannot be shifted forward because the' ...
                    ' frame-difference stage movement at ' ...
                    num2str(moveTime, '%.3f') ' seconds would have a' ...
                    ' predicted time of ' num2str(predictedTime, '%.3f') ...
                    ' seconds (an error of ' num2str(timeDiff, '%.3f') ...
                    ' seconds) whereas the interval between its media' ...
                    ' time and the previous media time is only ' ...
                    num2str(mediaDiff, '%.3f') ' seconds and,' ...
                    ' therefore, smaller than the error from shifting']);
                
            % Shift everything forward using the spare 0 media time location.
            elseif ~isempty(spareZeroTimeLocation)
                mediaTimes = [0 mediaTimes];
                locations = [spareZeroTimeLocation; locations];
                movesI(end + 1,:) = [0 0];
                timeOff = prevPeakI / fps - mediaTimes(i - 1);
                if verbose
                    peaksI = [1 peaksI];
                    otsuThrs = [gOtsuThr otsuThrs];
                    smallThrs = [gSmallThr smallThrs];
                    timeOffs(1:(end + 1)) = [0 timeOffs(1:end)];
                    timeDiffs =  movesI(2:(i - 2),1)' / fps - ...
                        mediaTimes(2:(i - 2));
                    timeOffs(3:(i - 1)) = mediaTimes(3:(i - 1)) + timeDiffs;
                end
                
                % Redo the match.
                i = i - 1;
                
                % Warn about the time shift.
                warning('findStageMovement:TimeShiftForward', ...
                    ['Shifting the media times forward relative to the ' ...
                    'frame-difference stage movements (using a spare ' ...
                    'location at media time 0:0:0.000) in an attempt ' ...
                    'to realign them']);
                
            % Shift everything forward by assuming a missing 0 media time
            % location and swallowing earlier frames into the the first
            % stage movement.
            else
                frame_statuses(1:movesI(2,1)) = true;
                movesI(1:(i - 2),:) = movesI(2:(i - 1),:);
                movesI(1,1) = 0;
                timeOff = prevPeakI / fps - mediaTimes(i - 2);
                if verbose
                    peaksI(1:(i - 2)) = peaksI(2:(i - 1));
                    otsuThrs(1:(i - 2)) = otsuThrs(2:(i - 1));
                    smallThrs(1:(i - 2)) = smallThrs(2:(i - 1));
                    timeDiffs =  movesI(1:(i - 3),1)' / fps - ...
                        mediaTimes(1:(i - 3));
                    timeOffs(2:(i - 2)) = mediaTimes(2:(i - 2)) + timeDiffs;
                end
                
                % Redo the match.
                i = i - 2;
                
                % Warn about the time shift.
                warning('findStageMovement:TimeShiftForward', ...
                    ['Shifting the media times forward relative to the ' ...
                    'frame-difference stage movements (by swallowing ' ...
                    'earlier frames into the first stage movement) in ' ...
                    'an attempt to realign them']);
            end
            
        % Shift the media times backward.
        else
            
            % Would a time shift align the media times with the
            % frame-difference stage movements?
            for j = 3:(i - 1)
                
                % Compute the error from the predicted time.
                offset =  movesI(j - 1,1) / fps - mediaTimes(j);
                predictedTime = mediaTimes(j + 1) + offset;
                moveTime =  movesI(j,1) / fps;
                timeDiff = abs(predictedTime - moveTime);
                
                % Compute the interval between the media times.
                mediaDiff = mediaTimes(j + 1) - mediaTimes(j);
                
                % Is the error in the predicted time greater than the
                % interval between media times?
                if timeDiff > mediaDiff
                    isShiftable = false;
                    break;
                end
            end
            
            % Time cannot be shifted due to misalignment between the media
            % times and frame-difference stage movements.
            if ~isShiftable
                warning('findStageMovement:TimeShiftAlignment', ...
                    ['Time cannot be shifted backward because the' ...
                    ' frame-difference stage movement at ' ...
                    num2str(moveTime, '%.3f') ' seconds would have a' ...
                    ' predicted time of ' num2str(predictedTime, '%.3f') ...
                    ' seconds (an error of ' num2str(timeDiff, '%.3f') ...
                    ' seconds) whereas the interval between its media' ...
                    ' time and the previous one is only ' ...
                    num2str(mediaDiff, '%.3f') ' seconds and,' ...
                    ' therefore, smaller than the error from shifting']);
                
            % Shift everything backward.
            else
                mediaTimes(1) = [];
                locations(1,:) = [];
                movesI(end,:) = [];
                timeOff = prevPeakI / fps - mediaTimes(i - 1);
                if verbose
                    peaksI(end) = [];
                    otsuThrs(end) = [];
                    smallThrs(end) = [];
                    timeOffs(1) = [];
                    timeDiffs =  movesI(1:(i - 2),1)' / fps - ...
                        mediaTimes(1:(i - 2));
                    timeOffs(1:(i - 1)) = [mediaTimes(1), ...
                        mediaTimes(2:(i - 1)) + timeDiffs];
                end
                
                % Redo the match.
                i = i - 1;
                
                % Warn about the time shift.
                warning('findStageMovement:TimeShiftBackward', ...
                    ['Shifting the media times backward relative to ' ...
                    'the frame-difference stage movements in an ' ...
                    'attempt to realign them']);
            end
        end
        
        % Record the shift and continue.
        if isShiftable
            isShifted = true;
            continue;
            
        % We cannot realign (shift) the stage movements and media times.
        else
            
            % Compute the stage movement sizes.
            movesI(i:end,:) = [];
            moveSizes = zeros(size(movesI, 1),1);
            for j = 2:(size(movesI, 1) - 1)
                moveDiffs = frameDiffs(movesI(j,1):movesI(j,2));
                moveSizes(j) = sum(moveDiffs(~isnan(moveDiffs)));
            end
            
            % Compute the statistics for stage movement sizes.
            meanMoveSize = mean(moveSizes(2:end));
            stdMoveSize = std(moveSizes(2:end));
            smallMoveThr = meanMoveSize - 2.5 * stdMoveSize;
            largeMoveThr = meanMoveSize + 2.5 * stdMoveSize;
            
            % Are any of the stage movements considerably small or large?
            for j = 2:(size(movesI, 1) - 1)
                
                % Is the stage movement small?
                if moveSizes(j) < smallMoveThr
                    
                    % Report the warning.
                    warning('findStageMovement:ShortMove', ...
                        ['Stage movement ' num2str(j) ...
                        ' at media time ' num2str(mediaTimes(j), '%.3f') ...
                        ' seconds (frame ' ...
                        num2str(round(mediaTimes(j) * fps)) ...
                        '), spanning from ' ...
                        num2str((movesI(j,1) - 1) / fps, '%.3f') ...
                        ' seconds (frame ' num2str(movesI(j,1) - 1) ...
                        ') to ' num2str((movesI(j,2) - 1) / fps, '%.3f') ...
                        ' seconds (frame ' ...
                        num2str(movesI(j,2) - 1) '), is considerably small']);
                    
                % Is the stage movement large?
                elseif moveSizes(j) > largeMoveThr
                    
                    % Report the warning.
                    warning('findStageMovement:LongMove', ...
                        ['Stage movement ' num2str(j) ...
                        ' at media time ' num2str(mediaTimes(j), '%.3f') ...
                        ' seconds (frame ' ...
                        num2str(round(mediaTimes(j) * fps)) ...
                        '), spanning from ' ...
                        num2str((movesI(j,1) - 1) / fps, '%.3f') ...
                        ' seconds (frame ' num2str(movesI(j,1) - 1) ...
                        ') to ' num2str((movesI(j,2) - 1) / fps, '%.3f') ...
                        ' seconds (frame ' ...
                        num2str(movesI(j,2) - 1) '), is considerably large']);
                end
            end
                    
            % Construct the report.
            id = 'findStageMovement:NoShift';
            msg = ['We cannot find a matching peak nor shift the time ' ...
                'for stage movement ' num2str(i) ' at media time ' ...
                num2str(mediaTimes(i), '%.3f') ' seconds (frame ' ...
                num2str(round(mediaTimes(i) * fps)) ')'];
        
            % Report the error.
            if ~verbose
                error(id, msg);
                
            % Report the warning.
            else
                warning(id, msg);
                
                % Finish.
                isMediaTimesMatched = false;
                peaksI(i:end) = [];
                otsuThrs(i:end) = [];
                smallThrs((i - 1):end) = [];
                timeOffs(i:end) = [];
                break;
            end
        end
    end
        
    % Find a temporary back end for this stage movement.
    % Note: this peak may serve as its own temporary back end.
    startI = max(peakI - maxMoveFrames, prevPeakEndI);
    [minDiff j] = min(fliplr(frameDiffs(startI:peakI)));
    peakBackEndI = peakI - j + 1; % we flipped to choose the last min
    j = peakI - 1;
    
    % If the temporary back end's frame difference is small, try to push
    % the back end forwards (closer to the stage movement).
    if minDiff <= prevSmallThr
        while j > startI
            if frameDiffs(j) <= prevSmallThr
                peakBackEndI = j;
                break;
            end
            j = j - 1;
        end
        
    % If the temporary back end's frame difference is large, try to push
    % the back end backwards (further from the stage movement).
    elseif minDiff >= min(otsuThr, gOtsuThr) || ...
            (minDiff > gSmallThr && peakBackEndI > startI && ...
            all(isnan(frameDiffs(startI:(peakBackEndI - 1)))))
        peakBackEndI = startI;
    end
    
    % Compute a threshold for stage movement.
    smallDiffs = frameDiffs(prevPeakEndI:peakBackEndI);
    smallThr = nanmean(smallDiffs) + 3 * nanstd(smallDiffs);
    if isnan(smallThr)
        smallThr = prevSmallThr;
    end
    if (verbose)
        smallThrs(i - 1) = smallThr;
    end
    
    % Find the front end for the previous stage movement.
    j = prevPeakI;
    while j < peakI && (isnan(frameDiffs(j)) || ...
            frameDiffs(j) > smallThr) && (isnan(frameDiffs(j + 1)) || ...
            frameDiffs(j + 1) > smallThr)
        j = j + 1;
    end
    movesI(i - 1,2) = j - 1;
    prevPeakEndI = j - 1;
    
    % Mark the previous stage movement.
    if movesI(i - 1,1) < 1
        frame_statuses(1:movesI(i - 1,2)) = true;
    else
        frame_statuses(movesI(i - 1,1):movesI(i - 1,2)) = true;
    end
    
    % Find the back end for this stage movement.
    j = peakI;
    while j > prevPeakEndI && (isnan(frameDiffs(j)) || ...
            frameDiffs(j) > smallThr)
        j = j - 1;
    end
    movesI(i, 1) = j + 1;
    
    % Is the non-movement frame-differences threshold too large?
    if smallThr <= otsuThr && (isnan(gOtsuThr) || smallThr <= gOtsuThr)
        prevOtsuThr = otsuThr;
        prevSmallThr = smallThr;
    else
        warning('findStageMovement:LargeNonMovementThreshold', ...
            ['The non-movement window between stage movement ' ...
            num2str(i - 1) ' and stage movement ' num2str(i) ...
            ', from ' num2str((movesI(i - 1,2) - 1) / fps, '%.3f') ...
            ' seconds (frame ' num2str(movesI(i - 1,2) - 1) ...
            ') to ' num2str((movesI(i,1) - 1) / fps, '%.3f') ...
            ' seconds (frame ' num2str(movesI(i,1) - 1) '),' ...
            ' contains considerably large frame-difference variance']);
    end
    
    % Compute the media time offset.
    timeOff = peakTime - mediaTimes(i);
    
    % We reached the end.
    endI = peakI + maxMoveFrames;
    if endI >= length(frameDiffs)
        peakFrontEndI = length(frameDiffs);
        
    % Find a temporary front end for this stage movement.
    else
        [minDiff j] = min(frameDiffs((peakI + 1):endI));
        peakFrontEndI = peakI + j;
        
        % If the temporary front end's frame difference is large, try to
        % push the front end forwards (further from the stage movement).
        if minDiff >= min(otsuThr, gOtsuThr) || ...
                (minDiff > max(smallThr, gSmallThr) && ...
                peakFrontEndI < endI && ...
                all(isnan(frameDiffs((peakFrontEndI + 1):endI))))
            peakFrontEndI = endI;
        end
    end
    
    % Try to push the temporary front end backwards (closer to the stage
    % movement).
    j = peakI + 1;
    while j < peakFrontEndI
        if frameDiffs(j) <= smallThr
            peakFrontEndI = j;
            break;
        end
        j = j + 1;
    end
    
    % Advance.
    prevPeakI = peakI;
    prevPeakEndI = peakFrontEndI;
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


end

function helper__adjustThreshold(gOtsuThr,gSmallThr,gSmallDiffs,frameDiffs)

%--------------------------------------------------------------------------
% Does the Otsu threshold separate the 99% of the small frame differences
% from the large ones?
if isempty(gSmallDiffs) || gOtsuThr < gSmallThr
    warning('findStageMovements:NoGlobalOtsuThreshold', ...
        ['Using the Otsu method, as a whole, the frame differences ' ...
        'don''t appear to contain any distinguishably large peaks ' ...
        '(corresponding to stage movements). Trying half of the ' ...
        'maximum frame difference instead.']);
    
    % Try half the maximum frame difference as a threshold to distinguish large peaks.
    gOtsuThr    = .5;
    gSmallDiffs = frameDiffs(frameDiffs < gOtsuThr);
    gSmallThr   = median(gSmallDiffs) + 3 * std(gSmallDiffs);
    
    % Does a threshold at half the maximum frame difference separate the
    % 99% of the small frame differences from the large ones?
    if isempty(gSmallDiffs) || gOtsuThr < gSmallThr
        warning('findStageMovements:NoGlobalThresholds', ...
            ['Cannot find a global threshold to distinguish the large ' ...
            'frame-difference peaks.']);
        gOtsuThr  = NaN;
        gSmallThr = NaN;
    end
    
end
%--------------------------------------------------------------------------




end

function helper__plotStuffs()

%{

if verbose
    
    % Open up a big figure.
    isGUI = true;
    if isempty(guiHandle)
        figure('OuterPosition', [50 50 1280 960]);
        guiHandle = axes;
        isGUI = false;
    end
    hold on;
    
    % Plot the video frame differences. Then plot the offset stage movement
    % times and their otsu thresholds on the same axes.
    [ax h1 h2] = plotyy(guiHandle, 0:(length(frameDiffs) - 1), ...
        frameDiffs, timeOffs, otsuThrs);
    set(ax(2), 'XAxisLocation', 'top');
    set(h1, 'Color', 'r');
    set(h2, 'Color', 'b', 'LineStyle', 'none', 'Marker', '.');
    
    % Setup the axes colors.
    set(ax(1), 'XColor', 'r', 'YColor', 'r');
    set(ax(2), 'XColor', 'k', 'YColor', 'k');
    
    % Setup the axes numbering.
    linkaxes(ax, 'y');
    xlim(ax(1), [0 (length(frameDiffs) - 1)]);
    xlim(ax(2), [0 ((length(frameDiffs) - 1) / fps)]);
    set(ax(2), 'YTick', []);
    set(zoom(guiHandle), 'Motion', 'horizontal', 'Enable', 'on');
    
    % Setup the axes labels.
    xlabel(ax(1), 'Frame');
    ylabel(ax(1), 'Subsequent Frame-Difference Variance');
    xlabel(ax(2), 'Time (seconds)');
    if ~isGUI % underscores confuse TeX
        title(ax(2), strrep(logFile, '_', '\_'));
    end
    
    % Setup the legends.
    legends1{1} = 'Variance of Subsequent Frame Differences';
    legends2 = [];
    if ~isempty(timeOffs)
        legends2{end + 1} = ...
            'Offset Movement Times at Otsu Threshold Height';
    end
    
    % Hold on.
    hold(ax(1), 'on');
    hold(ax(2), 'on');
    
    % Plot the offset stage movement times and their otsu thresholds.
    if verbose > 1
        text(timeOffs, otsuThrs, num2str((1:length(timeOffs))'), ...
            'Color', 'b', 'Parent', ax(2));
    end
    
    % Pretty plot the stage movements.
    % Note: earlier, we shifted the frame differences over one to align
    % them with the video frames. Due to differentiation and this shift,
    % the large differences at the end of stage movements represent a drop
    % in variance and, therefore, a non-movement frame. Visually, it's much
    % easier to recognize correct stage movement detection by aligning the
    % start and end of large differences with frame status. Therefore, we
    % extend the end of detected stage movements by one in order to achieve
    % visual alignment of large frame differences and stage movement intervals.
    plotFrames = frames;
    if movesI(1,2) > 0
        plotFrames(movesI(1,2) + 1) = true;
    end
    plotFrames(movesI(2:(end - 1),2) + 1) = true;
    plot(ax(1), 0:(length(plotFrames) - 1), plotFrames, 'k');
    if ~isempty(plotFrames)
        legends1{end + 1} = ...
            'Stage Movements (Shifted for Variance Alignment)';
    end
    
    % Plot the small thresholds.
    plotThrs(1:length(frames)) = NaN;
    startI = 1;
    if movesI(1,2) > 0
        startI = movesI(1,2) + 2;
    end
    plotThrs(startI:(movesI(2,1) - 1)) = smallThrs(1);
    for i = 2:(size(movesI, 1) - 1)
        plotThrs((movesI(i,2) + 2):(movesI(i + 1,1)) - 1) = smallThrs(i);
    end
    plot(ax(1), 0:(length(plotThrs) - 1), plotThrs, 'c');
    if ~isempty(plotThrs)
        legends1{end + 1} = 'Non-movement Thresholds';
    end
    
    % Plot the matched video frame difference peaks.
    brown = [.5 .3 .1];
    plot(ax(1), peaksI - 1, frameDiffs(peaksI), '.', 'Color', brown);
    if ~isempty(peaksI)
        legends1{end + 1} = 'Movement Peaks';
    end
    if verbose > 1
        matchedPeaksStr = num2str((1:length(peaksI))');
        text(peaksI - 1, frameDiffs(peaksI), matchedPeaksStr, ...
            'Color', brown, 'Parent', ax(1));
    end
    
    % Plot the unmatched video frame difference peaks.
    yellow = [1 .8 .2];
    plot(ax(1), endPeaksI - 1, frameDiffs(endPeaksI), '.', ...
        'Color', yellow);
    if ~isempty(endPeaksI)
        legends1{end + 1} = 'Unmatched Peaks';
    end
       
    % Plot the matched stage movement times and their small thresholds.
    matchedMediaTimes = mediaTimes(1:length(peaksI));
    matchedSmallThrs = [gSmallThr smallThrs(1:(length(peaksI) - 1))];
    plot(ax(2), matchedMediaTimes, matchedSmallThrs, 'g.');
    if ~isempty(matchedMediaTimes)
        legends2{end + 1} = ...
            'Movement Times at Non-Movement Threshold Height';
    end
    if verbose > 1
        text(matchedMediaTimes, matchedSmallThrs, matchedPeaksStr, ...
            'Color', 'g', 'Parent',ax(2));
    end
    
    % Plot the unmatched stage movement times.
    unmatchedMediaTimes = ...
        mediaTimes((length(peaksI) + 1):length(mediaTimes));
    unmatchedMediaTimesY = zeros(length(unmatchedMediaTimes), 1);
    plot(ax(2), unmatchedMediaTimes, unmatchedMediaTimesY, 'm.');
    if ~isempty(unmatchedMediaTimes)
        legends2{end + 1} = 'Unmatched Stage Movements';
    end
    if verbose > 1
        text(unmatchedMediaTimes, unmatchedMediaTimesY, ...
            num2str(((length(peaksI) + 1):length(mediaTimes))'), ...
            'Color', 'm', 'Parent', ax(2));
    end
    
    % Show the legend.
    legends1((end + 1):(end + length(legends2))) = legends2;
    legend(ax(1), legends1, 'Location', 'NorthEast');
    
    % Hold off.
    hold(ax(1), 'off');
    hold(ax(2), 'off');
    
    % Report the unmatched stage movements.
    if ~isempty(unmatchedMediaTimes)
        error('findStageMovement:UnmatchedStageMovements', ...
            ['There are ' num2str(length(unmatchedMediaTimes)) ...
            ' stage movements that couldn''t be matched to ' ...
            'large, appropriately-timed frame differences']);
    end
end

%}

end