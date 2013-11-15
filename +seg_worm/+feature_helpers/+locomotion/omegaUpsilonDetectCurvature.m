function [omegaFrames, upsilonFrames] = omegaUpsilonDetectCurvature(angleArray, stage_moving_flag)
% omegaUpsilonDetectCurvature  This function determines which frames of the input
% angle array contain part of an omega turn.  Omega turns have three
% characteristics: 
%
%
%   OLD NAME: 
%
%
%   1) They start with high curvature in the first third of the worm.
%   2) They continue with high curvature at the midbody.
%   3) They end with high curvature at the tail.
%
% Upsilon bends are similar to omega bends but have less pronounced
% curvature at each stage of the turn.
%
%   [omegaFrames, upsilonFrames] =
%       seg_worm.feature_helpers.locomotion.omegaUpsilonDetectCurvature(angleArray, stageFlag)
%
% Input:
%   angleChange   - matrix of skeleton tangent angle changes
%   stageFlag     - logical array indicating stage movement
%
%
% Output:
%   omegaFrames   - frames that contain part of an omega turn are indicated
%                   by +/- 1, other frames are 0. If the worm side is
%                   included, the sign indicates whether the bend is ventral
%                   (+1) or dorsal (-1)
%   upsilonFrames - frames that contain part of an upsilon turn are
%                   indicated by +/- 1, other frames are 0. If the worm
%                   side is included, the sign indicates whether the bend
%                   is ventral (+1) or dorsal (-1)
%

%Old code for validation ..
%[omegaFrames2, upsilonFrames2] = helper__omegaUpsilonDetectCurvature(angleArray, stage_moving_flag);

[numSegments, numFrames] = size(angleArray);

% initialize output arrays
f.omegaFrames   = zeros(numFrames, 1);
f.upsilonFrames = zeros(numFrames, 1);

FIRST_THIRD_INDICES  = 1:round(numSegments * (1/3));
SECOND_THIRD_INDICES = (FIRST_THIRD_INDICES(end)+1):round(numSegments*(2/3));
LAST_THIRD_INDICES   = (SECOND_THIRD_INDICES(end)+1):numSegments;

% divide angle array into three parts and take the mean of each
headAngle = nanmean(angleArray(FIRST_THIRD_INDICES, :));
bodyAngle = nanmean(angleArray(SECOND_THIRD_INDICES, :));
tailAngle = nanmean(angleArray(LAST_THIRD_INDICES, :));


%Possible early end ...
%--------------------------------------------------------------------------
n_head = sum(~isnan(headAngle));
n_body = sum(~isnan(bodyAngle));
n_tail = sum(~isnan(tailAngle));

%only proceed if there are at least two non-NaN value in each angle vector
if n_head < 2 || n_body < 2 || n_tail < 2
   return 
end

%Interpolation
%--------------------------------------------------------------------------

%Get long NaN stretches ...
n = isnan(bodyAngle);
%This little bit finds runs of NaN values that are 10 samples or more
%0 -> A
%1 -> B
[long_nan_start_I, long_nan_end_I] = regexp( char(n+'A'), 'B{10,}', 'start', 'end' );

% interpolate arrays over NaN values (where there were stage
% movements, touching, or some other segmentation problem)
% ***This is of course only an approximate solution to the problem of
% not segmenting coiled shapes***
a.headAngle = h__interp_NaN(headAngle);
a.bodyAngle = h__interp_NaN(bodyAngle);
a.tailAngle = h__interp_NaN(tailAngle);

% return long NaN stretches back to NaN
for kk = 1:length(long_nan_start_I)
    a.bodyAngle(long_nan_start_I(kk):long_nan_end_I(kk)) = NaN;
end

a.stage_moving_flag = stage_moving_flag;

%--------------------------------------------------------------------------

%This doesn't match was is written in the supplemental material ...
%Am I working off of old code??????
c = struct(...
    'headAngleStartConst',{20 -20 15 -15}, ...
    'tailAngleStartConst',{30 30  30  30}, ...
    'headAngleEndConst',  {40 40  30 30},  ...
    'tailAngleEndConst',  {20 -20 15 -15}, ...
    'bodyAngleConst'   ,  {20 -20 15 -15});

is_upsilon  = [false false true true];
sign_values = [1    -1     1    -1];

for iEntry = 1:4
    s = h__getStuffs(a,c(iEntry));
    f = h__populateFrames(a,s,f,is_upsilon(iEntry),sign_values(iEntry));
end

omegaFrames   = f.omegaFrames;
upsilonFrames = f.upsilonFrames;

end

function f = h__populateFrames(a,s,f,get_upsilon_flag,sign_value)
    
    %Algorithm:
    %-----------------------------------------------------------
    %- For the middle angle range, ensure one frame is valid and that
    %  the frame proceeding the start and following the end are valid
    %- Find start indices and end indices that bound this range
    %- For upsilons, exclude if they overlap with an omega bend ...
    
    %JAH NOTE: This type of searching is inefficient in Matlab since 
    %the data is already sorted. It could be improved ...
    
    for iMid = 1:length(s.midStarts)
        cur_mid_start_I = s.midStarts(iMid);
        cur_mid_end_I   = s.midEnds(find(s.midEnds > cur_mid_start_I, 1));
        
        if ~isempty(cur_mid_end_I) && ...
            ~all(a.stage_moving_flag(cur_mid_start_I:cur_mid_end_I)) && ...
            s.startCond(cur_mid_start_I - 1) && ...
            s.endCond(cur_mid_end_I + 1)
 
            cur_start_I = s.startInds(find(s.startInds < cur_mid_start_I,          1, 'last'));
            cur_end_I   = s.endInds(find(s.endInds     > cur_mid_end_I, 1));

            if get_upsilon_flag
                %Don't populate upsilon if the data spans an omega
                if ~any(abs(f.omegaFrames(cur_start_I:cur_end_I)))
                    f.upsilonFrames(cur_start_I:cur_end_I) = sign_value;
                end
            else
                f.omegaFrames(cur_start_I:cur_end_I) = sign_value;
            end
        end
    end
end


function s = h__getStuffs(a,c)

    is_positive = c.headAngleStartConst > 0;

    if is_positive
        fh = @gt;
    else
        fh = @lt;
    end
    
    %NOTE:
    %start: when the head exceeds its angle but the tail does not
    %end  : when the tail exceeds its angle but the head does not
    
    s.startCond = fh(a.headAngle, c.headAngleStartConst) & abs(a.tailAngle) < c.tailAngleStartConst;
    s.startInds = find(diff(s.startCond) == 1) + 1; %add 1 for shift due to diff
    
    %NOTE: This is NaN check is a bit suspicious, as it implies that the
    %head and tail are parsed, but the body is not. The original code puts
    %NaN back in for long gaps in the body angle, so it is possible that
    %the body angle is NaN but the others are not.
    s.midCond   = fh(a.bodyAngle, c.bodyAngleConst) | isnan(a.bodyAngle);
    s.midStarts = find(diff(s.midCond) == 1) + 1; %add 1 for shift due to diff
    s.midEnds   = find(diff(s.midCond) == -1);
    
    s.endCond   = fh(a.tailAngle, c.tailAngleEndConst) & abs(a.headAngle) < c.headAngleEndConst;
    s.endInds   = find(diff(s.endCond) == -1);

end


%This is old code ...
%------------------------------------------------------------

%{
function [omegaFrames, upsilonFrames] = helper__omegaUpsilonDetectCurvature(angleArray, stageFlag)
% OMEGAUPSILONDETECTDV This function determines which frames of the input
% angle array contain part of an omega turn.  Omega turns have three
% characteristics: they start with high curvature in the first third of the
% worm, they continue with high curvature at the midbody, and end with high
% curvature at the tail.
%
%
% Upsilon bends are similar to omega bends but have less pronounced
% curvature at each stage of the turn.
%
% Input:
%   angleChange   - matrix of skeleton tangent angle changes
%   stageFlag     - logical array indicating stage movement
%   tailToHeadDirectionBlock
% Output:
%   omegaFrames   - frames that contain part of an omega turn are indicated
%                   by +/- 1, other frames are 0. If the worm side is
%                   included, the sign indicates whether the bend is ventral
%                   (+1) or dorsal (-1)
%   upsilonFrames - frames that contain part of an upsilon turn are
%                   indicated by +/- 1, other frames are 0. If the worm
%                   side is included, the sign indicates whether the bend
%                   is ventral (+1) or dorsal (-1)
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copysright notices and other proprietary 
% notices on any copies of the Software.

[numSegments, numFrames] = size(angleArray);

% initialize output arrays
omegaFrames   = zeros(numFrames, 1);
upsilonFrames = zeros(numFrames, 1);

FIRST_THIRD_INDICES  = 1:round(numSegments * (1/3));
SECOND_THIRD_INDICES = (FIRST_THIRD_INDICES(end)+1):round(numSegments*(2/3));
LAST_THIRD_INDICES   = (SECOND_THIRD_INDICES(end)+1):numSegments;

% divide angle array into three parts and take the mean of each
headAngle = nanmean(angleArray(FIRST_THIRD_INDICES, :));
bodyAngle = nanmean(angleArray(SECOND_THIRD_INDICES, :));
tailAngle = nanmean(angleArray(LAST_THIRD_INDICES, :));

n_head = sum(~isnan(headAngle));
n_body = sum(~isnan(bodyAngle));
n_tail = sum(~isnan(tailAngle));

%only proceed if there are at least two non-NaN value in each angle vector
if n_head < 2 || n_body < 2 || n_tail < 2
   return 
end


    
    % interpolate arrays over NaN values (where there were stage
    % movements, touching, or some other segmentation problem)
    % ***This is of course only an approximate solution to the problem of
    % not segmenting coiled shapes***
    headAngle(isnan(headAngle)) = interp1(find(~isnan(headAngle)),...
        headAngle(~isnan(headAngle)), find(isnan(headAngle)),'linear', 'extrap');
    tailAngle(isnan(tailAngle)) = interp1(find(~isnan(tailAngle)),...
        tailAngle(~isnan(tailAngle)), find(isnan(tailAngle)),'linear', 'extrap');
    
    % Interpolate body angle
    % get long NaN stretches
    n = isnan(bodyAngle); % reshape
    % save start and end indices for the stretches
    [start1, end1] = regexp( char(n+'A'), 'B{10,}', 'start', 'end' );
    
    bodyAngle(isnan(bodyAngle)) = interp1(find(~isnan(bodyAngle)),...
        bodyAngle(~isnan(bodyAngle)), find(isnan(bodyAngle)),'linear', 'extrap');
    
    % return long NaN stretches back to NaN
    if ~isempty(start1) && ~isempty(end1)
        for kk=1:size(start1,2)
            bodyAngle(start1(kk):end1(kk)) = NaN;
        end
    end
%==========================================================================
%==========================================================================
    
    
    
    
    
    
    
    %********************** Find Positive Omegas **************************
    
    headAngleStartConst = 20;
    headAngleEndConst   = 40;
    tailAngleStartConst = 30;
    tailAngleEndConst   = 20;
    bodyAngleConst      = 20;
    
    %find frames that satisfy conditions for omega bend start
    startCond = headAngle > headAngleStartConst & abs(tailAngle) < tailAngleStartConst;
    startInds = find(diff(startCond) == 1) + 1; %add 1 for shift due to diff
    
    %find frames that satisfy conditions for middle of omega bend (these
    %could be coils or stage motions, so include check for NaNs
    midCond   = bodyAngle > bodyAngleConst | isnan(bodyAngle);
    
    midStarts = find(diff(midCond) == 1) + 1; %add 1 for shift due to diff
    midEnds   = find(diff(midCond) == -1);
    
    %find frames that satisfy conditions for omega bend end
    endCond = tailAngle > tailAngleEndConst & abs(headAngle) < headAngleEndConst;
    endInds = find(diff(endCond) == -1);
    
    for j = 1:length(midStarts)
        % find the next end index that is greater than the current startInd
        possibleEnd = find(midEnds > midStarts(j), 1);
        
        % stage motions are allowed during turns, but there must be at
        % least one valid non-motion frame.
        if all(stageFlag(midStarts(j):midEnds(possibleEnd)))
            continue;
        end
                
        % check that frames before and after the possible omega turn are
        % valid start and end frames respectively
        if ~isempty(possibleEnd)
            startCheck1 = midStarts(j)-1;
            startCheck2 = midStarts(j)-1;
            
            endCheck1 = midEnds(possibleEnd)+1;
            endCheck2 = midEnds(possibleEnd)+1;
            if startCheck1 > 0 && endCheck2 <= length(endCond)
                if any(startCond(startCheck1:startCheck2)) && any(endCond(endCheck1:endCheck2))
                    % we have a positive omega turn.  Now find actual start and end
                    % points.
                    currentStart = find(startInds < midStarts(j), 1, 'last');
                    currentEnd = find(endInds > midEnds(possibleEnd), 1);
                    
                    omegaFrames(startInds(currentStart):endInds(currentEnd)) = 1;
                end
            end
        end
    end
    
    
    
    
    
    
    
    
    
    
    %********************** Find Negative Omegas **************************
    headAngleStartConst = -20;
    headAngleEndConst   = 40;
    tailAngleStartConst = 30;
    tailAngleEndConst   = -20;
    bodyAngleConst      = -20;
    
    %find frames that satisfy conditions for omega bend start
    startCond = headAngle < headAngleStartConst & abs(tailAngle) < tailAngleStartConst;
    startInds = find(diff(startCond) == 1) + 1; %add 1 for shift due to diff
    
    %find frames that satisfy conditions for middle of omega bend
    midCond = bodyAngle < bodyAngleConst | isnan(bodyAngle);
    midStarts = find(diff(midCond) == 1) + 1; %add 1 for shift due to diff
    midEnds = find(diff(midCond) == -1);
    
    %find frames that satisfy conditions for omega bend end
    endCond = tailAngle < tailAngleEndConst & abs(headAngle) < headAngleEndConst;
    endInds = find(diff(endCond) == -1);
    
    for j = 1:length(midStarts)
        % find the next end index that is greater than the current startInd
        possibleEnd = find(midEnds > midStarts(j), 1);
        
        % stage motions are allowed during turns, but there must be at
        % least one valid non-motion frame.
        if all(stageFlag(midStarts(j):midEnds(possibleEnd)))
            continue;
        end
                
        % check that frames before and after the possible omega turn are
        % valid start and end frames respectively
        if ~isempty(possibleEnd)
            startCheck1 = midStarts(j)-1;
            startCheck2 = midStarts(j)-1;
            
            endCheck1 = midEnds(possibleEnd)+1;
            endCheck2 = midEnds(possibleEnd)+1;
            if startCheck1 > 0 && endCheck2 <= length(endCond)
                if any(startCond(startCheck1:startCheck2)) && any(endCond(endCheck1:endCheck2))
                    % we have a positive omega turn.  Now find actual start and end
                    % points.
                    currentStart = find(startInds < midStarts(j), 1, 'last');
                    currentEnd = find(endInds > midEnds(possibleEnd), 1);
                    
                    omegaFrames(startInds(currentStart):endInds(currentEnd)) = -1;
                end
            end
        end
    end
    
    
    
    
    
    
    
    
    %********************** Find Positive Upsilons ************************
    headAngleStartConst = 15;
    headAngleEndConst   = 30;
    tailAngleStartConst = 30;
    tailAngleEndConst   = 15;
    bodyAngleConst      = 15;
    
    %find frames that satisfy conditions for upsilon bend start
    startCond = headAngle > headAngleStartConst & abs(tailAngle) < tailAngleStartConst;
    startInds = find(diff(startCond) == 1) + 1; %add 1 for shift due to diff
    
    %find frames that satisfy conditions for middle of upsilon bend
    midCond = bodyAngle > bodyAngleConst | isnan(bodyAngle);
    midStarts = find(diff(midCond) == 1) + 1; %add 1 for shift due to diff
    midEnds = find(diff(midCond) == -1);
    
    %find frames that satisfy conditions for upsilon bend end
    endCond = abs(headAngle) < headAngleEndConst & tailAngle > tailAngleEndConst;
    endInds = find(diff(endCond) == -1);
    
    for j = 1:length(midStarts)
        % find the next end index that is greater than the current startInd
        possibleEnd = find(midEnds > midStarts(j), 1);
        
        % stage motions are allowed during turns, but there must be at
        % least one valid non-motion frame.
        if all(stageFlag(midStarts(j):midEnds(possibleEnd)))
            continue;
        end
                
        % check that frames before and after the possible upsilon turn are
        % valid start and end frames respectively
        if ~isempty(possibleEnd)
            if startCond(midStarts(j) - 1) && endCond(midEnds(possibleEnd) + 1)
                % we have a positive upsilon turn.  Now find actual start and end
                % points.
                currentStart = find(startInds < midStarts(j), 1, 'last');
                currentEnd = find(endInds > midEnds(possibleEnd), 1);
                
                % Here we will check if an upsilon bend has an omega bend
                % inside. If yes we will not count it and will go to the
                % next iteration.
                if ~any(abs(omegaFrames( ...
                        startInds(currentStart):endInds(currentEnd) )))
                    upsilonFrames( ...
                         startInds(currentStart):endInds(currentEnd) ) = 1;
                end
            end
        end
    end
    
    %********************** Find Negative Upsilons ************************
    headAngleStartConst = -15;
    headAngleEndConst   = 30;
    tailAngleStartConst = 30;
    tailAngleEndConst   = -15;
    bodyAngleConst      = -15;
    
    %find frames that satisfy conditions for upsilon bend start
    startCond = headAngle < headAngleStartConst & abs(tailAngle) < tailAngleStartConst;
    startInds = find(diff(startCond) == 1) + 1; %add 1 for shift due to diff
    
    %find frames that satisfy conditions for middle of upsilon bend
    midCond   = bodyAngle < bodyAngleConst | isnan(bodyAngle);
    midStarts = find(diff(midCond) == 1) + 1; %add 1 for shift due to diff
    midEnds   = find(diff(midCond) == -1);
    
    %find frames that satisfy conditions for upsilon bend end
    endCond = abs(headAngle) < headAngleEndConst & tailAngle < tailAngleEndConst;
    endInds = find(diff(endCond) == -1);
    
    for j = 1:length(midStarts)
        % find the next end index that is greater than the current startInd
        possibleEnd = find(midEnds > midStarts(j), 1);
        
        % stage motions are allowed during turns, but there must be at
        % least one valid non-motion frame.
        if all(stageFlag(midStarts(j):midEnds(possibleEnd)))
            continue;
        end
        
        % check that frames before and after the possible upsilon turn are
        % valid start and end frames respectively
        if ~isempty(possibleEnd)
            if startCond(midStarts(j) - 1) && endCond(midEnds(possibleEnd) + 1)
                % we have a positive upsilon turn.  Now find actual start and end
                % points.
                currentStart = find(startInds < midStarts(j), 1, 'last');
                currentEnd   = find(endInds > midEnds(possibleEnd), 1);
                
                % Here we will check if an upsilon bend has an omega bend
                % inside. If yes we will not count it and will go to the
                % next iteration.
                if ~any(abs(omegaFrames(startInds(currentStart):endInds(currentEnd) )))
                    upsilonFrames(startInds(currentStart):endInds(currentEnd) ) = -1;
                end
            end
        end
    end
end

function fixed_x = h__interp_NaN(x)

fixed_x  = x;
nan_mask = isnan(x);

fixed_x(nan_mask) = interp1(find(~nan_mask),x(~nan_mask), find(nan_mask),'linear', 'extrap');

end
%}
